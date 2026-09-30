#!/usr/bin/env python3
"""Write one surface field from a NetCDF file as a WPS intermediate file.

This writes a single 2-D field on a regular latitude-longitude grid, under the
name init_atmosphere looks for, to a file named PREFIX:YYYY-MM-DD_hh:

    ./write_field_intermediate.py era5_lsm.nc 2000-01-01_00:00:00 \\
        --var lsm --field LANDSEA --prefix LSM --match ERA5:2000-01-01_00

produces LSM:2000-01-01_00. The typical use is to supply a field that is
missing from the ungribbed meteorological data, such as the LANDSEA mask, by
appending the file to it:

    cat LSM:2000-01-01_00 >> ERA5:2000-01-01_00

init_atmosphere reads LANDSEA at the grid indices of each field it
interpolates, so the mask must be on exactly the same grid, in the same order,
as the rest of the file. With --match, the output is oriented to the grid of an
existing intermediate file (ungrib often writes latitude north to south), and
the script stops if the two grids differ. Without it, latitude runs south to
north and longitude from 0 to 360, as in write_sst_intermediate.py.

Masked, NaN and fill values become -1.e30, the missing value. The file is read
back after writing and summarized, as a check on the format.
"""

import argparse
import struct
import sys

import numpy as np
from netCDF4 import Dataset

from write_sst_intermediate import (TIME_RE, MISSING, write_field, read_back,
                                    read_2d, regular_step)


def read_grid(path):
    """Return (nx, ny, startlat, startlon, dlat, dlon) of the first field in an
    intermediate file, which must be on a lat-lon grid."""
    def rec(f):
        n, = struct.unpack('>i', f.read(4))
        payload = f.read(n)
        f.read(4)
        return payload

    with open(path, 'rb') as f:
        version, = struct.unpack('>i', rec(f))
        if version != 5:
            sys.exit('{}: expected format version 5, got {}.'.format(path, version))
        h = rec(f)
        nx, ny, iproj = struct.unpack('>3i', h[144:156])
        if iproj != 0:
            sys.exit('{}: the first field is not on a lat-lon grid (iproj = {}); '
                     'only lat-lon grids are supported.'.format(path, iproj))
        g = rec(f)
        startlat, startlon, dlat, dlon = struct.unpack('>4f', g[8:24])
    return nx, ny, startlat, startlon, dlat, dlon


def match_grid(path, data, lat, lon):
    """Reorder data (on lat south to north, lon in [0, 360)) to the grid of the
    intermediate file at path, or stop if the grids differ."""
    nx, ny, startlat, startlon, dlat, dlon = read_grid(path)
    if (ny, nx) != data.shape:
        sys.exit('The grid is {} x {} (lat x lon), but {} is {} x {}.'
                 .format(data.shape[0], data.shape[1], path, ny, nx))
    tol = 1.e-3 * max(abs(dlat), abs(dlon))
    if dlat < 0:
        lat, data = lat[::-1], data[::-1, :]
    # Start the longitudes at startlon, wrapping around the globe.
    offset = (lon - startlon + 180.) % 360. - 180.
    shift = int(np.argmin(np.abs(offset)))
    if abs(offset[shift]) > tol:
        sys.exit('No longitude of the grid is at {:g}, the start longitude of {}: the '
                 'nearest is {:g}.'.format(startlon, path, lon[shift]))
    lon, data = np.roll(lon, -shift), np.roll(data, -shift, axis=1)
    lon = startlon + np.mod(lon - startlon + tol, 360.) - tol

    ours = (lat[0], lon[0], lat[1] - lat[0], lon[1] - lon[0])
    theirs = (startlat, startlon, dlat, dlon)
    if not np.allclose(ours, theirs, rtol=0., atol=tol):
        sys.exit('The grid (startlat, startlon, dlat, dlon) = ({:g}, {:g}, {:g}, {:g}) '
                 'does not match {}: ({:g}, {:g}, {:g}, {:g}).'
                 .format(*ours, path, *theirs))
    print('Matched the grid of {}.'.format(path))
    return data, lat, lon


def find_coord(nc, given, candidates):
    if given:
        return given
    for name in candidates:
        if name in nc.variables:
            return name
    sys.exit('None of {} is in the file; name the coordinate with --lat or --lon.'
             .format(', '.join(candidates)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file', help='NetCDF file with the field')
    p.add_argument('time', help='valid time, e.g. 2000-01-01_00:00:00, written in the '
                                'header and used in the output file name')
    p.add_argument('--var', required=True, help='name of the NetCDF variable')
    p.add_argument('--field', required=True,
                   help='name of the field in the intermediate file, e.g. LANDSEA')
    p.add_argument('--prefix', required=True,
                   help='output prefix: the file is written as PREFIX:YYYY-MM-DD_hh')
    p.add_argument('--match', metavar='INTERMEDIATE_FILE',
                   help='existing intermediate file whose grid and ordering the '
                        'output must match, e.g. ERA5:2000-01-01_00')
    p.add_argument('--units', help='units to write (default: the units attribute)')
    p.add_argument('--desc', help='description to write (default: the long_name attribute)')
    p.add_argument('--lat', help="latitude coordinate (default 'lat' or 'latitude')")
    p.add_argument('--lon', help="longitude coordinate (default 'lon' or 'longitude')")
    p.add_argument('--time-index', type=int, default=0,
                   help='index along the time axis, if there is one (default 0)')
    args = p.parse_args()

    if not TIME_RE.match(args.time):
        sys.exit('The time must have the form YYYY-MM-DD_hh:mm:ss: {} does not.'
                 .format(args.time))
    if len(args.field) > 9:
        sys.exit('The field name can be at most 9 characters: {} is not.'.format(args.field))

    nc = Dataset(args.file)
    lat_name = find_coord(nc, args.lat, ('lat', 'latitude'))
    lon_name = find_coord(nc, args.lon, ('lon', 'longitude'))
    lat = np.asarray(nc.variables[lat_name][:], dtype='f8')
    lon = np.asarray(nc.variables[lon_name][:], dtype='f8')
    if lat.ndim != 1 or lon.ndim != 1:
        sys.exit('Latitude and longitude must be 1-D coordinates of a regular grid.')
    data, units = read_2d(nc, args.var, lat_name, lon_name, args.time_index)
    # The header has fixed-width slots: 25 characters for units, 46 for the
    # description.
    desc = (args.desc or getattr(nc.variables[args.var], 'long_name', args.field))[:46]
    units = (args.units or units)[:25]
    nc.close()

    # South to north.
    if lat[0] > lat[-1]:
        lat, data = lat[::-1], data[::-1, :]
    # Longitude in [0, 360), increasing, so the south-west corner and a positive
    # spacing describe the grid.
    lon = np.mod(lon, 360.)
    order = np.argsort(lon, kind='stable')
    lon, data = lon[order], data[:, order]

    regular_step(lat, 'Latitude')
    dlon = regular_step(lon, 'Longitude')
    if abs(dlon * lon.size - 360.) > 0.5 * dlon:
        sys.exit('The grid covers {:g} degrees of longitude; init_atmosphere wraps '
                 'the field around the globe, so it must cover 360.'.format(dlon * lon.size))

    if args.match:
        data, lat, lon = match_grid(args.match, data, lat, lon)

    out = '{}:{}'.format(args.prefix, args.time[:13])
    slab = np.where(np.isfinite(data), data, MISSING)
    with open(out, 'wb') as f:
        write_field(f, args.time, args.field, units, desc, slab,
                    lat[0], lon[0], lat[1] - lat[0], lon[1] - lon[0])

    print('Wrote {} ({} x {} grid, dlat = {:g}, dlon = {:g}). Read back:'
          .format(out, lat.size, lon.size, lat[1] - lat[0], lon[1] - lon[0]))
    read_back(out)


if __name__ == '__main__':
    main()
