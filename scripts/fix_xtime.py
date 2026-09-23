#!/usr/bin/env python3
"""Set the 'xtime' variable of an MPAS file.

MPAS seeks on 'xtime' when it reads a stream, so a file whose 'xtime' is blank,
NUL-padded or truncated fails with

    ERROR: Invalid DateTime string (invalid date substring)

and the string printed after the message is empty, because there is nothing in
it to print. This writes a well-formed, blank-padded timestamp in place.

    ./fix_xtime.py x1.2562.init.nc 2000-01-01_00:00:00

The timestamp must match config_start_time exactly: the initial read seeks on
an exact time, not the nearest one.
"""

import argparse
import re
import sys

import numpy as np
from netCDF4 import Dataset

TIME_RE = re.compile(r'^\d{4,}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2}$')


def show(label, raw):
    """Print a char record the way MPAS sees it, with invisible bytes named."""
    text = ''.join(c.decode('ascii', 'replace') if isinstance(c, bytes) else str(c)
                   for c in np.asarray(raw).ravel())
    # A NUL, tab or carriage return is what usually breaks the file, and each is
    # invisible on its own, so show every one of them as '@' and name its code.
    shown = ''.join(c if 32 <= ord(c) <= 126 else '@' for c in text)
    print('  {:8s} [{}]'.format(label, shown))
    odd = sorted({ord(c) for c in text if not (32 <= ord(c) <= 126)})
    if odd:
        print('  {:8s} the @ marks are ASCII {}'
              .format('', ', '.join(str(c) for c in odd)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file')
    p.add_argument('time', help='timestamp, e.g. 2000-01-01_00:00:00')
    p.add_argument('--record', type=int, default=None,
                   help='Time record to set (default: every record)')
    p.add_argument('--dry-run', action='store_true',
                   help='report what is in the file and stop')
    args = p.parse_args()

    if not TIME_RE.match(args.time):
        sys.exit('The timestamp must have the form YYYY-MM-DD_hh:mm:ss, with every '
                 'field padded to its full width: {} does not.'.format(args.time))

    nc = Dataset(args.file, 'r' if args.dry_run else 'r+')
    # Without this, netCDF4 masks NUL bytes away and the padding that broke the
    # file becomes invisible in the report below.
    nc.set_auto_mask(False)

    if 'xtime' in nc.variables:
        xtime = nc.variables['xtime']
        strlen = xtime.shape[-1]
    else:
        # The variable is absent, so create it alongside the Time dimension.
        if args.dry_run:
            sys.exit("'xtime' is not in the file; rerun without --dry-run to create it.")
        if 'Time' not in nc.dimensions:
            nc.createDimension('Time', None)
        if 'StrLen' not in nc.dimensions:
            nc.createDimension('StrLen', 64)
        strlen = len(nc.dimensions['StrLen'])
        xtime = nc.createVariable('xtime', 'S1', ('Time', 'StrLen'))
        print("Created 'xtime' with StrLen = {}".format(strlen))

    if len(args.time) > strlen:
        sys.exit('The timestamp is {} characters but StrLen is {}.'
                 .format(len(args.time), strlen))

    nrecords = xtime.shape[0]
    records = range(nrecords) if args.record is None else [args.record]

    print('{}: {} Time record(s), StrLen = {}'.format(args.file, nrecords, strlen))
    for i in records:
        if i < nrecords:
            show('was', xtime[i, :])

    if args.dry_run:
        nc.close()
        return

    # Blank padding, not NUL padding: trim() in MPAS treats a NUL as a character.
    padded = np.array(list(args.time.ljust(strlen)), dtype='S1')
    for i in records:
        xtime[i, :] = padded
        show('now', xtime[i, :])

    nc.close()
    print('Done. Check with: ncdump -v xtime {}'.format(args.file))


if __name__ == '__main__':
    main()
