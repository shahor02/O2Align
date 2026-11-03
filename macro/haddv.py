#!/usr/bin/env python3
"""hadd with concatenation of std::vector<o2::dataformats::TFIDInfo>.

Run inside an O2 environment with PyROOT and hadd on PATH.
Examples:
  ./haddv.py merged.root input1.root input2.root
  ./haddv.py -f merged.root @inputs.txt
  ./haddv.py -a merged.root more/*.root

Supports the common hadd flags listed by --help, rather than every hadd flag.
Vectors (including those in directories) are concatenated in source order,
without sorting or deduplication. Only the latest key cycle is used.
Local output only; input ROOT URLs are supported by ROOT.
"""

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


def arguments(argv):
    # hadd encodes compression in the force flag, e.g. -f505 or -fk206.
    compression = []
    normalized = []
    for arg in argv:
        if re.fullmatch(r'-f(?:\d+|k\d*|f)', arg):
            compression.append(arg)
            normalized.append('-f')
        else:
            normalized.append(arg)
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    for flag, help_text in [('-f', 'overwrite output (also -fN, -fk[N], -ff)'),
                            ('-a', 'append to existing output'),
                            ('-k', 'skip unreadable input files'),
                            ('-O', 're-optimize tree baskets'),
                            ('-T', 'do not merge trees')]:
        p.add_argument(flag, action='store_true', help=help_text)
    for flag in ('-v', '-n', '-cachesize'):
        p.add_argument(flag, metavar='VALUE')
    p.add_argument('target')
    p.add_argument('sources', nargs='+', help='ROOT files or @file lists (one path per line)')
    args = p.parse_args(normalized)
    args.compression = compression[-1:] or (['-f'] if args.f else [])
    return args


def expand_sources(sources, stack=()):
    for source in sources:
        if not source.startswith('@'):
            yield source
            continue
        path = Path(source[1:]).resolve()
        if path in stack:
            raise ValueError(f'Recursive input list: {path}')
        with path.open() as stream:
            entries = [line.strip() for line in stream if line.strip() and not line.lstrip().startswith('#')]
        yield from expand_sources(entries, (*stack, path))


def open_root(ROOT, name, mode='READ'):
    f = ROOT.TFile.Open(str(name), mode)
    if not f or f.IsZombie():
        if f:
            f.Close()
        raise OSError(f'Cannot open ROOT file: {name}')
    return f


def vectors(ROOT, directory, prefix=''):
    seen = set()
    for key in directory.GetListOfKeys():
        name = key.GetName()
        if name in seen:
            continue
        seen.add(name)
        cls = ROOT.TClass.GetClass(key.GetClassName())
        path = prefix + name
        if cls and cls.InheritsFrom('TDirectory'):
            yield from vectors(ROOT, directory.Get(name), path + '/')
        elif key.GetClassName().replace('std::', '').replace(' ', '') == 'vector<o2::dataformats::TFIDInfo>':
            yield path, directory.Get(name)


def main(argv=None):
    args = arguments(sys.argv[1:] if argv is None else argv)
    if '://' in args.target:
        raise ValueError('The output must be a local file')
    target = Path(args.target).absolute()
    if target.exists() and not (args.f or args.a):
        raise ValueError(f'Output exists (use -f to overwrite): {target}')
    sources = list(expand_sources(args.sources))
    if not sources:
        raise ValueError('No input files')
    for source in sources:
        if '://' not in source and (Path(source).resolve() == target.resolve() or
                                   (target.exists() and Path(source).exists() and os.path.samefile(source, target))):
            raise ValueError(f'Output is also an input: {source}')
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)
    vector_type = ROOT.std.vector('o2::dataformats::TFIDInfo')
    merged = {}
    valid = []
    if args.a and target.exists():
        sources.insert(0, str(target))
    for source in sources:
        try:
            f = open_root(ROOT, source)
        except (OSError, ReferenceError) as error:
            if not args.k or source == str(target):
                raise
            print(f'Skipping {source}: {error}', file=sys.stderr)
            continue
        try:
            for path, vector in vectors(ROOT, f):
                if path not in merged:
                    merged[path] = vector_type()
                merged[path].insert(merged[path].end(), vector.begin(), vector.end())
            valid.append(source)
        finally:
            f.Close()
    if not valid:
        raise ValueError('No readable inputs')
    # Keep the old destination intact until both merging steps have succeeded.
    with tempfile.TemporaryDirectory(prefix='.hadd_tfidinfo-', dir=target.parent) as work:
        temporary = Path(work) / 'merged.root'
        options = list(args.compression)
        for flag in ('O', 'T'):
            if getattr(args, flag):
                options.append('-' + flag)
        for flag in ('v', 'n', 'cachesize'):
            value = getattr(args, flag)
            if value is not None:
                options.extend(['-' + flag, value])
        # A list avoids command-line length limits with large sets of inputs.
        listing = Path(work) / 'inputs.txt'
        listing.write_text(''.join(source + '\n' for source in valid))
        subprocess.run([shutil.which('hadd') or 'hadd', *options, str(temporary), '@' + str(listing)], check=True)
        output = open_root(ROOT, temporary, 'UPDATE')
        try:
            for path, vector in merged.items():
                parent, _, name = path.rpartition('/')
                directory = output.GetDirectory(parent) if parent else output
                if not directory:
                    raise RuntimeError(f'Missing output directory for {path}')
                if directory.WriteObject(vector, name, 'Overwrite') <= 0:
                    raise OSError(f'Failed to write {path}')
                print(f'{path}: {vector.size()} TFIDInfo records')
        finally:
            output.Close()
        if output.TestBit(ROOT.TFile.kWriteError):
            raise OSError('ROOT reported an output write error')
        if args.f or args.a:
            os.replace(temporary, target)
        else:
            # Fail safely if another process created the target meanwhile.
            os.link(temporary, target)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f'hadd_tfidinfo: {error}', file=sys.stderr)
        sys.exit(1)
