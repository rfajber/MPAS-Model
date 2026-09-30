#!/usr/bin/env python3
"""Write SST and sea ice from a NetCDF file as a WPS intermediate file.

init_atmosphere reads the sea surface temperature and sea-ice fraction for the
initial conditions from a separate intermediate file when config_input_sst =
true (see src/core_init_atmosphere/README.separate_sst_seaice.md). This writes
that file from SST and sea-ice data on a regular latitude-longitude grid:

    ./write_sst_intermediate.py oisst.nc 2020-01-01_00:00:00 \\
        --sst-var sst --ice-var ice

produces SST:2020-01-01_00, with fields named SST (K) and SEAICE (fraction), in
the version-5 format that mpas_init_atm_read_met.F reads.

The input is tidied so that init_atmosphere interprets it correctly:
  * SST in degrees Celsius is converted to kelvin. init_atmosphere turns any
    water cell below 271 K into sea ice, so an SST left in Celsius would cover
    the ocean in ice.
  * Sea ice in percent is converted to a fraction and clipped to [0, 1].
  * Masked, NaN and fill values become -1.e30, the missing value, so that
    init_atmosphere searches for the nearest valid point instead of
    interpolating from them.
  * Latitude is flipped to run south to north, and longitude is shifted to
    [0, 360) and sorted, since the format describes the grid by its south-west
    corner and spacing.

With --mask-ice-with-sst, the points where the SST is missing are taken as land
and the sea ice is masked with them: it is set to missing there, and the mask is
written as LANDSEA, which init_atmosphere uses to exclude land points when it
interpolates SEAICE. Sea ice is then interpolated only from points that have a
valid SST, even if the sea-ice data set has its own, different land mask.

The grid must be regular and global in longitude. The file is read back after
writing and summarized, as a check on the format.
"""

import argparse
import re
import struct
import sys

import numpy as np
from netCDF4 import Dataset

TIME_RE = re.compile(r'^\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}$')
MISSING = -1.e30            # msgval in mpas_init_atm_surface.F
SURFACE_LEVEL = 200100.     # the intermediate-format code for a surface field
EARTH_RADIUS_KM = 6367.47   # written for completeness; unused for lat-lon grids


def _rec(f, payload):
    """Write one big-endian Fortran unformatted sequential record."""
    n = len(payload)
    f.write(struct.pack('>i', n))
    f.write(payload)
    f.write(struct.pack('>i', n))


def write_field(f, hdate, field, units, desc, data, startlat, startlon, dlat, dlon):
    """Append one field in version-5 format. data is (nlat, nlon), with
    latitude increasing from startlat and longitude increasing from startlon."""
    ny, nx = data.shape
    _rec(f, struct.pack('>i', 5))
    _rec(f, hdate.ljust(24).encode()
            + struct.pack('>f', 0.0)                   # forecast hour
            + 'write_sst_intermediate'.ljust(32).encode()
            + field.ljust(9).encode()
            + units.ljust(25).encode()
            + desc.ljust(46).encode()
            + struct.pack('>f', SURFACE_LEVEL)
            + struct.pack('>3i', nx, ny, 0))           # iproj 0 = lat-lon
    _rec(f, b'SWCORNER'
            + struct.pack('>5f', startlat, startlon, dlat, dlon, EARTH_RADIUS_KM))
    _rec(f, struct.pack('>i', 0))                      # is_wind_grid_rel = .false.
    # C order on (nlat, nlon) puts longitude fastest, which is Fortran slab(nx, ny).
    _rec(f, np.ascontiguousarray(data, dtype='>f4').tobytes())


def read_back(path):
    """Read the file the way mpas_init_atm_read_met.F does and summarize it."""
    def rec(f):
        head = f.read(4)
        if not head:
            return None
        n, = struct.unpack('>i', head)
        payload = f.read(n)
        tail, = struct.unpack('>i', f.read(4))
        if tail != n:
            raise ValueError('record markers disagree ({} vs {})'.format(n, tail))
        return payload

    with open(path, 'rb') as f:
        while True:
            r = rec(f)
            if r is None:
                break
            version, = struct.unpack('>i', r)
            if version != 5:
                raise ValueError('expected format version 5, got {}'.format(version))
            h = rec(f)
            hdate = h[0:24].decode().strip()
            field = h[60:69].decode().strip()
            units = h[69:94].decode().strip()
            xlvl, = struct.unpack('>f', h[140:144])
            nx, ny, iproj = struct.unpack('>3i', h[144:156])
            g = rec(f)
            startloc = g[0:8].decode()
            startlat, startlon, dlat, dlon, _ = struct.unpack('>5f', g[8:28])
            rec(f)
            slab = np.frombuffer(rec(f), dtype='>f4').reshape(ny, nx)
            valid = slab[slab > MISSING / 10]
            print('  {:8s} {:10s} {} nx={} ny={} iproj={} level={:g} {} '
                  'lat0={:g} lon0={:g} dlat={:g} dlon={:g}'
                  .format(field, units, hdate, nx, ny, iproj, xlvl, startloc,
                          startlat, startlon, dlat, dlon))
            if valid.size:
                print('  {:8s} valid {:.4g} to {:.4g}, {:.1f}% missing'
                      .format('', valid.min(), valid.max(),
                              100. * (1 - valid.size / slab.size)))
            else:
                print('  {:8s} every point is missing'.format(''))


def regular_step(values, name):
    """Return the spacing of a 1-D coordinate, or stop if it is irregular."""
    step = np.diff(values)
    if not np.allclose(step, step[0], rtol=1e-4, atol=1e-6):
        sys.exit('{} is not regularly spaced (steps from {:g} to {:g}); only regular '
                 'lat-lon grids are supported.'.format(name, step.min(), step.max()))
    return float(step[0])


def read_2d(nc, var, lat_name, lon_name, time_index):
    """Return a variable as a (nlat, nlon) float64 array with missing values as NaN."""
    v = nc.variables[var]
    dims = list(v.dimensions)
    if lat_name not in dims or lon_name not in dims:
        sys.exit("'{}' has dimensions {}; expected '{}' and '{}' among them."
                 .format(var, dims, lat_name, lon_name))
    index = []
    for d in dims:
        if d in (lat_name, lon_name):
            index.append(slice(None))
        elif v.shape[dims.index(d)] == 1:
            index.append(0)
        else:
            index.append(time_index)   # the time (or other) axis
    data = np.ma.filled(np.ma.masked_invalid(v[tuple(index)]).astype('f8'), np.nan)
    kept = [d for d in dims if d in (lat_name, lon_name)]
    if kept == [lon_name, lat_name]:
        data = data.T
    return data, getattr(v, 'units', '')


def to_kelvin(sst, units):
    u = units.strip().lower()
    celsius = u in ('c', 'degc', 'deg_c', 'celsius', 'degrees_c', 'degree_c',
                    'degrees_celsius', 'deg c', 'degrees c') or 'celsius' in u
    kelvin = u in ('k', 'kelvin', 'degk', 'deg_k', 'degrees_k')
    if not celsius and not kelvin:
        # No usable units attribute: decide from the values, and say so.
        celsius = np.nanmax(sst) < 100.
        print("SST units '{}' not recognized; values suggest {}."
              .format(units, 'Celsius' if celsius else 'kelvin'))
    if celsius:
        print('Converting SST from Celsius to kelvin.')
        sst = sst + 273.15
    return sst


def to_fraction(ice, units):
    u = units.strip().lower()
    if u in ('%', 'percent', 'percentage'):
        print('Converting sea ice from percent to a fraction.')
        ice = ice / 100.
    elif np.nanmax(ice) > 1.5:
        print("Sea ice reaches {:g} with units '{}'; treating it as percent."
              .format(np.nanmax(ice), units))
        ice = ice / 100.
    return np.clip(ice, 0., 1.)


def mask_ice_with_sst(data):
    """Take the points where the SST is missing as land, and apply that mask to
    the sea ice in both of the ways init_atmosphere recognizes: the ice is set to
    missing there, and the mask is written as LANDSEA, which interp_sfc_to_MPAS
    uses to exclude source land points when it interpolates SEAICE."""
    land = ~np.isfinite(data['SST'])
    ice = data['SEAICE']
    # Some SST products leave the SST missing under sea ice. Those points would
    # be treated as land and their ice dropped, so report them.
    lost = land & np.isfinite(ice) & (ice > 0.)
    if lost.any():
        print('Warning: {} points have sea ice but no SST, and will be treated as land. '
              'If your SST product masks the SST under ice, do not use '
              '--mask-ice-with-sst.'.format(int(lost.sum())))
    print('Masking sea ice with the SST mask: {} of {} points are land.'
          .format(int(land.sum()), land.size))
    data['SEAICE'] = np.where(land, np.nan, ice)
    data['LANDSEA'] = land.astype('f8')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file', help='NetCDF file with SST (and optionally sea ice)')
    p.add_argument('time', help='valid time, e.g. 2020-01-01_00:00:00; must equal '
                                'config_start_time in namelist.init_atmosphere')
    p.add_argument('--sst-var', required=True, help='name of the SST variable')
    p.add_argument('--ice-var', help='name of the sea-ice variable. Without it no '
                                     'SEAICE is written and init_atmosphere keeps '
                                     'the sea ice from the met file')
    p.add_argument('--landsea-var', help='optional land mask (1 = land, 0 = water), '
                                         'written as LANDSEA. init_atmosphere uses it '
                                         'to exclude land points from SEAICE only')
    p.add_argument('--mask-ice-with-sst', action='store_true',
                   help='treat every point where the SST is missing as land: set the '
                        'sea ice to missing there and write the mask as LANDSEA, so '
                        'that init_atmosphere interpolates SEAICE only from points '
                        'with a valid SST. Needs --ice-var; replaces --landsea-var')
    p.add_argument('--lat', default='lat', help="latitude coordinate (default 'lat')")
    p.add_argument('--lon', default='lon', help="longitude coordinate (default 'lon')")
    p.add_argument('--time-index', type=int, default=0,
                   help='index along the time axis to write (default 0)')
    p.add_argument('--prefix', default='SST',
                   help="output prefix, i.e. config_sfc_prefix (default 'SST')")
    args = p.parse_args()

    if not TIME_RE.match(args.time):
        sys.exit('The time must have the form YYYY-MM-DD_hh:mm:ss: {} does not.'
                 .format(args.time))
    if args.mask_ice_with_sst and not args.ice_var:
        sys.exit('--mask-ice-with-sst needs --ice-var: there is no sea ice to mask.')
    if args.mask_ice_with_sst and args.landsea_var:
        sys.exit('Use either --landsea-var or --mask-ice-with-sst, not both: each '
                 'supplies the LANDSEA field.')

    nc = Dataset(args.file)
    lat = np.asarray(nc.variables[args.lat][:], dtype='f8')
    lon = np.asarray(nc.variables[args.lon][:], dtype='f8')
    if lat.ndim != 1 or lon.ndim != 1:
        sys.exit('Latitude and longitude must be 1-D coordinates of a regular grid.')

    fields = [('SST', args.sst_var), ('SEAICE', args.ice_var), ('LANDSEA', args.landsea_var)]
    data = {}
    for name, var in fields:
        if var:
            data[name], units = read_2d(nc, var, args.lat, args.lon, args.time_index)
            if name == 'SST':
                data[name] = to_kelvin(data[name], units)
            elif name == 'SEAICE':
                data[name] = to_fraction(data[name], units)
    nc.close()

    if args.mask_ice_with_sst:
        mask_ice_with_sst(data)

    # South to north.
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        data = {k: v[::-1, :] for k, v in data.items()}
    # Longitude in [0, 360), increasing, so the south-west corner and a positive
    # spacing describe the grid.
    lon = np.mod(lon, 360.)
    order = np.argsort(lon, kind='stable')
    lon = lon[order]
    data = {k: v[:, order] for k, v in data.items()}

    dlat = regular_step(lat, 'Latitude')
    dlon = regular_step(lon, 'Longitude')
    if abs(dlon * lon.size - 360.) > 0.5 * dlon:
        sys.exit('The grid covers {:g} degrees of longitude; init_atmosphere wraps '
                 'the field around the globe, so it must cover 360.'.format(dlon * lon.size))

    sst = data['SST']
    ocean = np.isfinite(sst)
    if ocean.any() and np.nanmin(sst) < 250.:
        print('Warning: SST as low as {:.1f} K. Check the units and fill values: water '
              'cells below 271 K become sea ice.'.format(np.nanmin(sst)))

    out = '{}:{}'.format(args.prefix, args.time[:13])
    desc = {'SST': 'Sea surface temperature', 'SEAICE': 'Sea-ice fraction',
            'LANDSEA': 'Land-sea mask (1=land)'}
    units = {'SST': 'K', 'SEAICE': 'fraction', 'LANDSEA': 'proprtn'}
    with open(out, 'wb') as f:
        for name, _ in fields:
            if name in data:
                slab = np.where(np.isfinite(data[name]), data[name], MISSING)
                write_field(f, args.time, name, units[name], desc[name], slab,
                            lat[0], lon[0], dlat, dlon)

    print('Wrote {} ({} x {} grid, dlat = {:g}, dlon = {:g}). Read back:'
          .format(out, lat.size, lon.size, dlat, dlon))
    read_back(out)


if __name__ == '__main__':
    main()
