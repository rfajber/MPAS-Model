#!/usr/bin/env python3
"""Compare a field in two WPS intermediate files.

Checks whether two intermediate files hold the same SST data, or the same data
for any other field named with --field:

    ./compare_intermediate.py SST:2000-01-01_00 OTHER:2000-01-01_00
    ./compare_intermediate.py ERA5C:2000-01-01_00 SST:2000-01-01_00 --field SST --field SEAICE

For each field and level it compares the grid and then the values: the points
missing in one file and not the other, and the largest and mean difference
between the points valid in both, with where the largest one is on a lat-lon
grid. Values within --tol of each other count as equal; the default, 0,
requires them to be identical. The dates in the headers are reported but not
compared, so a field written for a different time can still match.

Two files with the same values on grids in a different order (latitude north
to south in one and south to north in the other, say) are reported as having
different grids: the values are compared point by point, as init_atmosphere
would read them.

The exit status is 0 if every field compared is the same, and 1 otherwise.
"""

import argparse
import struct
import sys

import numpy as np

MISSING = -1.e30


def rec(f):
    """Return the payload of one Fortran unformatted record, or None at the end."""
    head = f.read(4)
    if not head:
        return None
    n, = struct.unpack('>i', head)
    payload = f.read(n)
    f.read(4)
    return payload


def read_fields(path):
    """Return {(field, level): dict} for the first occurrence of each field and
    level in an intermediate file, warning about later ones."""
    fields = {}
    with open(path, 'rb') as f:
        while True:
            version = rec(f)
            if version is None:
                break
            if struct.unpack('>i', version)[0] != 5:
                sys.exit('{}: expected format version 5, got {}.'
                         .format(path, struct.unpack('>i', version)[0]))
            h, grid, _, slab = rec(f), rec(f), rec(f), rec(f)
            name = h[60:69].decode().strip()
            level, = struct.unpack('>f', h[140:144])
            nx, ny, iproj = struct.unpack('>3i', h[144:156])
            key = (name, level)
            if key in fields:
                print('Warning: {} has more than one {} at level {:g}; using the first.'
                      .format(path, name, level))
                continue
            fields[key] = {
                'date': h[0:24].decode().strip(),
                'units': h[69:94].decode().strip(),
                'nx': nx, 'ny': ny, 'iproj': iproj,
                'grid': grid,
                'data': np.frombuffer(slab, dtype='>f4').reshape(ny, nx).astype('f8'),
            }
    return fields


def latlon_grid(field):
    """Return (startlat, startlon, dlat, dlon) for a lat-lon field, else None."""
    if field['iproj'] != 0:
        return None
    return struct.unpack('>4f', field['grid'][8:24])


def describe_grid(field):
    g = latlon_grid(field)
    if g is None:
        return '{} x {}, iproj {}'.format(field['ny'], field['nx'], field['iproj'])
    return ('{} x {}, startlat {:g}, startlon {:g}, dlat {:g}, dlon {:g}'
            .format(field['ny'], field['nx'], *g))


def compare(name, level, a, b, tol):
    """Print how field a differs from field b, and return True if they match."""
    print('{} at level {:g}:'.format(name, level))
    if a['date'] != b['date']:
        print('  dates differ ({} and {}); not counted as a difference'
              .format(a['date'], b['date']))
    if a['units'] != b['units']:
        print("  units differ: '{}' and '{}'".format(a['units'], b['units']))
    if (a['nx'], a['ny'], a['iproj'], a['grid']) != (b['nx'], b['ny'], b['iproj'], b['grid']):
        print('  DIFFERENT GRIDS, values not compared:')
        print('    file 1: ' + describe_grid(a))
        print('    file 2: ' + describe_grid(b))
        return False

    va, vb = a['data'] > MISSING / 10, b['data'] > MISSING / 10
    miss_1, miss_2, both = int((vb & ~va).sum()), int((va & ~vb).sum()), va & vb
    diff = np.abs(a['data'] - b['data'])
    diff[~both] = 0.
    worst = np.unravel_index(np.argmax(diff), diff.shape)
    max_diff = diff[worst]
    over = int((diff > tol).sum())

    print('  {} points, {} missing in both, {} missing only in file 1, {} only in file 2'
          .format(diff.size, int((~va & ~vb).sum()), miss_1, miss_2))
    if both.any():
        print('  over the {} points valid in both: max |difference| {:.6g}, mean {:.6g}, '
              '{} differ by more than {:g}'
              .format(int(both.sum()), max_diff, diff[both].mean(), over, tol))
        if max_diff > 0:
            j, i = worst
            where = 'row {}, column {}'.format(j + 1, i + 1)
            g = latlon_grid(a)
            if g is not None:
                where += ' (lat {:g}, lon {:g})'.format(g[0] + j * g[2], g[1] + i * g[3])
            print('  largest at {}: {:.9g} in file 1, {:.9g} in file 2'
                  .format(where, a['data'][worst], b['data'][worst]))

    same = miss_1 == 0 and miss_2 == 0 and over == 0 and a['units'] == b['units']
    if same and np.array_equal(a['data'], b['data']):
        print('  IDENTICAL')
    elif same:
        print('  SAME within {:g}'.format(tol))
    else:
        print('  DIFFERENT')
    return same


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file1', help='first intermediate file')
    p.add_argument('file2', help='second intermediate file')
    p.add_argument('--field', action='append',
                   help="field to compare; repeat for several (default 'SST')")
    p.add_argument('--tol', type=float, default=0.,
                   help='largest difference that counts as equal (default 0, exact)')
    args = p.parse_args()
    names = args.field or ['SST']

    a, b = read_fields(args.file1), read_fields(args.file2)
    all_same = True
    for name in names:
        levels = sorted({lev for n, lev in a if n == name} | {lev for n, lev in b if n == name})
        if not levels:
            print('{}: in neither file.'.format(name))
            all_same = False
        for level in levels:
            if (name, level) not in a or (name, level) not in b:
                print('{} at level {:g}: only in {}.'.format(
                    name, level, args.file1 if (name, level) in a else args.file2))
                all_same = False
            else:
                all_same = compare(name, level, a[(name, level)], b[(name, level)],
                                   args.tol) and all_same

    print('Result: {}'.format('same' if all_same else 'different'))
    sys.exit(0 if all_same else 1)


if __name__ == '__main__':
    main()
