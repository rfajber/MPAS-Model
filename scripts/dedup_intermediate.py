#!/usr/bin/env python3
"""Remove duplicate fields from a WPS intermediate file.

init_atmosphere stops with "Attempting to allocate already allocated variable
'maskslab'" when the meteorological file holds more than one LANDSEA field,
which happens when the same field is appended twice (for example with
write_field_intermediate.py and cat). This copies a file, keeping only the
first occurrence of each field at each level:

    ./dedup_intermediate.py ERA5F:2000-01-01_00 ERA5C:2000-01-01_00

Records are copied unchanged. Every version-5 field is five records (version,
header, grid, wind flag, data), so this works for any projection.

With --soilhgt-to-m, a SOILHGT field that holds geopotential (m2 s-2) rather
than height is divided by g and relabelled m. init_atmosphere does not convert
units, and would read geopotential as terrain about 9.8 times too high. The
field is taken to be geopotential if its units say so, or if it exceeds 9000,
which no terrain height does.

On a regional domain, check the result. Geopotential that is mislabelled m is
recognized only by its values, and over terrain lower than about 900 m it stays
below 9000 and is left unchanged. The script prints the units and maximum it
found: a maximum about 9.8 times the highest terrain in the domain means the
field is still geopotential.
"""

import argparse
import re
import struct

import numpy as np

G = 9.80665
MISSING = -1.e30
MAX_TERRAIN_M = 9000.   # above Everest: larger values can only be geopotential


def is_geopotential_units(units):
    """True for m**2 s**-2, m2 s-2, m^2/s^2 and the like."""
    return re.sub(r'[\s*^]', '', units.lower()) in ('m2s-2', 'm2/s2')


def rec(f):
    """Return one Fortran unformatted record, markers included, or None at the end."""
    head = f.read(4)
    if not head:
        return None
    n, = struct.unpack('>i', head)
    return head + f.read(n + 4)


def pack(payload):
    n = struct.pack('>i', len(payload))
    return n + payload + n


def soilhgt_to_m(recs):
    """Return the records of a SOILHGT field with the data in m, converting
    from geopotential if needed."""
    h = recs[1][4:-4]
    units = h[69:94].decode().strip()
    data = np.frombuffer(recs[4][4:-4], dtype='>f4').astype('f8')
    valid = data > MISSING / 10
    top = data[valid].max() if valid.any() else 0.
    by_units = is_geopotential_units(units)
    if top <= MAX_TERRAIN_M and not by_units:
        print("SOILHGT is in m (units '{}', maximum {:.1f}); left unchanged."
              .format(units, top))
        return recs
    reason = ("its units are '{}'".format(units) if by_units else
              'it reaches {:.1f}, above any terrain height'.format(top))
    data[valid] /= G
    print('SOILHGT taken as geopotential because {}; divided by g. It now '
          'reaches {:.1f} m.'.format(reason, top / G))
    h = h[:69] + 'm'.ljust(25).encode() + h[94:]
    return [recs[0], pack(h), recs[2], recs[3],
            pack(data.astype('>f4').tobytes())]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('src', help='intermediate file to read')
    p.add_argument('dst', help='intermediate file to write')
    p.add_argument('--soilhgt-to-m', action='store_true',
                   help='convert SOILHGT from geopotential (m2 s-2) to height (m) '
                        'if it is not in m already')
    args = p.parse_args()

    seen, kept, dropped = set(), 0, 0
    with open(args.src, 'rb') as f, open(args.dst, 'wb') as g:
        while True:
            recs = [rec(f) for _ in range(5)]   # version, header, grid, wind flag, data
            if recs[0] is None:
                break
            h = recs[1][4:-4]
            key = (h[60:69].decode().strip(), struct.unpack('>f', h[140:144])[0])
            if key in seen:
                dropped += 1
                print('dropping duplicate {} at level {:g}'.format(*key))
                continue
            seen.add(key)
            kept += 1
            if args.soilhgt_to_m and key[0] == 'SOILHGT':
                recs = soilhgt_to_m(recs)
            g.write(b''.join(recs))
    print('kept {} fields, dropped {}'.format(kept, dropped))


if __name__ == '__main__':
    main()
