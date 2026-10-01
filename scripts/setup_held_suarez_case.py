#!/usr/bin/env python3
"""Set up a case directory for the Held-Suarez (1994) benchmark on the 120-km mesh.

The case runs the dry dynamical core with the Held-Suarez forcing in place of
physics (see src/core_atmosphere/dynamics/README.held_suarez.md) on the uniform
x1.40962 mesh, starting at 2000-01-01_00:00:00 and running for 180 days:

    ./setup_held_suarez_case.py /path/to/x1.40962 \\
        --case-dir hs_120km --output-dir /scratch/hs_120km

The grid directory must hold x1.40962.grid.nc, and normally x1.40962.graph.info
and the x1.40962.graph.info.part.* files for the processor counts to be used.

Files that are used unchanged are linked into the case directory: the two
executables, the grid and partition files, streams.init_atmosphere and the
stream_list.atmosphere.* files. Files that are edited are copied from the
default_inputs directory of a built MPAS tree and then changed:

  namelist.init_atmosphere
      config_init_case = 2 (dry Jablonowski-Williamson baroclinic wave with a
      perturbation), the start time and config_nvertlevels. The JW case fixes
      its own 45-km model top and uniform level spacing, so config_ztop has no
      effect.
  namelist.atmosphere
      the start time, a 180-day run, config_dt = 720 and config_len_disp =
      120000 for the 120-km mesh, config_physics_suite = 'none' and
      config_held_suarez = .true.. The other config_hs_* options are left at
      the published benchmark values.
  streams.atmosphere
      daily history output, restarts every 30 days and no diagnostics stream,
      which carries mostly physics fields. With --output-dir, the history,
      diagnostics and restart files are written there instead of in the case
      directory.

Each option is replaced where it appears, and inserted into its record when it
is not in the defaults; a missing record or stream is an error, so that a
change to the Registry fails here rather than leaving a setting silently
unapplied.

Nothing is run. Afterwards, in the case directory, run init_atmosphere_model to
make x1.40962.init.nc and then atmosphere_model.
"""

import argparse
import glob
import os
import re
import shutil
import sys

MESH = 'x1.40962'
START_TIME = '2000-01-01_00:00:00'
RUN_DURATION = '180_00:00:00'
DT = '720.0'                     # 6 s per km of mesh spacing
LEN_DISP = '120000.0'            # m, the nominal mesh spacing
HISTORY_INTERVAL = '1_00:00:00'
RESTART_INTERVAL = '30_00:00:00'

INIT_NAMELIST_OPTIONS = {
    'nhyd_model': {
        'config_init_case': '2',
        'config_start_time': f"'{START_TIME}'",
    },
    'dimensions': {
        'config_nvertlevels': None,          # set from --nvertlevels
    },
}

ATM_NAMELIST_OPTIONS = {
    'nhyd_model': {
        'config_dt': DT,
        'config_start_time': f"'{START_TIME}'",
        'config_run_duration': f"'{RUN_DURATION}'",
        'config_len_disp': LEN_DISP,
    },
    'physics': {
        'config_physics_suite': "'none'",
    },
    'held_suarez': {
        'config_held_suarez': '.true.',
    },
}

# stream name -> {attribute: value}
ATM_STREAM_ATTRIBUTES = {
    'output': {'output_interval': HISTORY_INTERVAL},
    'restart': {'output_interval': RESTART_INTERVAL},
    'diagnostics': {'output_interval': 'none'},
}

# Streams whose files are written to --output-dir when it is given
OUTPUT_STREAMS = ['output', 'diagnostics', 'restart']

EXECUTABLES = ['init_atmosphere_model', 'atmosphere_model']
LINKED_INPUTS = [
    'streams.init_atmosphere',
    'stream_list.atmosphere.output',
    'stream_list.atmosphere.diagnostics',
    'stream_list.atmosphere.surface',
    'stream_list.atmosphere.diag_ugwp',
]
COPIED_INPUTS = [
    'namelist.init_atmosphere',
    'namelist.atmosphere',
    'streams.atmosphere',
]


def fail(msg):
    sys.exit('Error: ' + msg)


def set_namelist_options(text, options, fname):
    """Set each option in its record, inserting it before the record's closing
    '/' when the defaults leave it out."""
    for record, opts in options.items():
        m = re.search(r'^&' + re.escape(record) + r'\s*\n(.*?)^/', text,
                      flags=re.MULTILINE | re.DOTALL)
        if m is None:
            fail(f'record &{record} not found in {fname}')
        body = m.group(1)
        for name, value in opts.items():
            line = f'    {name} = {value}\n'
            pattern = r'^[ \t]*' + re.escape(name) + r'[ \t]*=.*\n'
            if re.search(pattern, body, flags=re.MULTILINE | re.IGNORECASE):
                body = re.sub(pattern, lambda _: line, body,
                              flags=re.MULTILINE | re.IGNORECASE)
            else:
                body += line
        text = text[:m.start(1)] + body + text[m.end(1):]
    return text


def set_stream_attribute(text, stream, attr, value, fname):
    """Set an attribute in the opening tag of a stream, keeping the layout the
    generated file already has."""
    m = re.search(r'<(?:immutable_)?stream\s+name="' + re.escape(stream) + r'"[^>]*>',
                  text)
    if m is None:
        fail(f'stream "{stream}" not found in {fname}')
    tag = m.group(0)
    attr_re = r'(\s' + re.escape(attr) + r'\s*=\s*)"[^"]*"'
    if re.search(attr_re, tag):
        tag = re.sub(attr_re, lambda a: f'{a.group(1)}"{value}"', tag)
    else:
        indent = re.search(r'\n(\s*)\S', tag)
        sep = '\n' + indent.group(1) if indent else ' '
        tag = re.sub(r'\s*/?>$', lambda e: f'{sep}{attr}="{value}"{e.group(0)}', tag)
    return text[:m.start()] + tag + text[m.end():]


def get_stream_attribute(text, stream, attr, fname):
    m = re.search(r'<(?:immutable_)?stream\s+name="' + re.escape(stream) + r'"[^>]*>',
                  text)
    if m is None:
        fail(f'stream "{stream}" not found in {fname}')
    a = re.search(r'\s' + re.escape(attr) + r'\s*=\s*"([^"]*)"', m.group(0))
    if a is None:
        fail(f'stream "{stream}" in {fname} has no {attr}')
    return a.group(1)


def is_same_link(src, dst):
    return os.path.islink(dst) and os.path.realpath(dst) == os.path.realpath(src)


def edit_atm_streams(output_dir):
    def edit(text, fname):
        for stream, attrs in ATM_STREAM_ATTRIBUTES.items():
            for attr, value in attrs.items():
                text = set_stream_attribute(text, stream, attr, value, fname)
        if output_dir:
            for stream in OUTPUT_STREAMS:
                template = get_stream_attribute(text, stream, 'filename_template', fname)
                template = os.path.join(output_dir, os.path.basename(template))
                text = set_stream_attribute(text, stream, 'filename_template',
                                            template, fname)
        return text
    return edit


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('grid_dir', help=f'directory holding {MESH}.grid.nc and the '
                   'graph.info partition files')
    p.add_argument('--case-dir', default=f'held_suarez_{MESH}',
                   help=f"case directory to create (default 'held_suarez_{MESH}')")
    p.add_argument('--output-dir',
                   help='directory the model writes its history, diagnostics and '
                   'restart files to (default: the case directory)')
    p.add_argument('--nvertlevels', type=int, default=26,
                   help='number of vertical levels (default 26)')
    p.add_argument('--mpas-root',
                   default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')),
                   help='built MPAS-Model tree providing the executables and '
                   'default_inputs (default: this repository)')
    p.add_argument('--overwrite', action='store_true',
                   help='replace files and links that already exist in the case directory')
    args = p.parse_args()

    if args.nvertlevels < 2:
        fail('--nvertlevels must be at least 2')

    mpas_root = os.path.abspath(args.mpas_root)
    grid_dir = os.path.abspath(args.grid_dir)
    case_dir = os.path.abspath(args.case_dir)
    output_dir = os.path.abspath(args.output_dir) if args.output_dir else None

    #
    # Collect the files to link and copy, and check that they all exist
    #
    links = {}     # name in case directory -> source
    copies = {}
    missing = []

    for f in EXECUTABLES:
        links[f] = os.path.join(mpas_root, f)
    for f in LINKED_INPUTS:
        links[f] = os.path.join(mpas_root, 'default_inputs', f)
    for f in COPIED_INPUTS:
        copies[f] = os.path.join(mpas_root, 'default_inputs', f)
    for src in list(links.values()) + list(copies.values()):
        if not os.path.isfile(src):
            missing.append(src)
    if missing:
        fail('these files were not found; has MPAS been built for both the '
             f'init_atmosphere and atmosphere cores in {mpas_root}?\n  '
             + '\n  '.join(missing))

    grid_file = os.path.join(grid_dir, f'{MESH}.grid.nc')
    if not os.path.isfile(grid_file):
        fail(f'{grid_file} not found')
    links[os.path.basename(grid_file)] = grid_file

    graph_files = sorted(glob.glob(os.path.join(grid_dir, f'{MESH}.graph.info*')))
    for g in graph_files:
        links[os.path.basename(g)] = g
    if not any('.part.' in g for g in graph_files):
        print(f'Warning: no {MESH}.graph.info.part.* files in {grid_dir}; '
              'only serial runs will be possible', file=sys.stderr)

    #
    # Refuse to touch an existing case directory's files unless asked to
    #
    if not args.overwrite:
        conflicts = [n for n, src in links.items()
                     if os.path.lexists(os.path.join(case_dir, n))
                     and not is_same_link(src, os.path.join(case_dir, n))]
        conflicts += [n for n in copies if os.path.lexists(os.path.join(case_dir, n))]
        if conflicts:
            fail(f'these files already exist in {case_dir}; use --overwrite to '
                 'replace them:\n  ' + '\n  '.join(conflicts))

    os.makedirs(case_dir, exist_ok=True)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    #
    # Link the files that are used unchanged
    #
    for name, src in links.items():
        dst = os.path.join(case_dir, name)
        if is_same_link(src, dst):
            continue
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(src, dst)
        print(f'Linked  {name}')

    #
    # Copy and edit the rest
    #
    init_options = {rec: dict(opts) for rec, opts in INIT_NAMELIST_OPTIONS.items()}
    init_options['dimensions']['config_nvertlevels'] = str(args.nvertlevels)

    edits = {
        'namelist.init_atmosphere':
            lambda t, f: set_namelist_options(t, init_options, f),
        'namelist.atmosphere':
            lambda t, f: set_namelist_options(t, ATM_NAMELIST_OPTIONS, f),
        'streams.atmosphere': edit_atm_streams(output_dir),
    }

    for name, src in copies.items():
        dst = os.path.join(case_dir, name)
        if os.path.lexists(dst):
            os.remove(dst)
        shutil.copyfile(src, dst)
        with open(dst) as f:
            text = f.read()
        text = edits[name](text, dst)
        with open(dst, 'w') as f:
            f.write(text)
        print(f'Copied  {name} (edited)')

    print()
    print(f'Held-Suarez case set up in {case_dir}')
    print(f'  mesh {MESH}, {args.nvertlevels} levels, start {START_TIME}, '
          f'run duration {RUN_DURATION}')
    print(f'  model output written to {output_dir or case_dir}')
    print('Next, in the case directory:')
    print('  1. run init_atmosphere_model to create '
          f'{MESH}.init.nc')
    print('  2. run atmosphere_model')


if __name__ == '__main__':
    main()
