#!/usr/bin/env python3
"""Copy selected GRID files, preserving their directory tree.

Use -d for comma-separated GRID folders or a text file with one folder per line,
and -l / --files for comma-separated file names. The alien:// prefix is optional.
Paths are replicated beneath -o (an existing directory; default: current directory).
Use -n to process multiple folders concurrently (default: 1), and -p to run an
executable in each folder after all its copies succeed. Existing requested files
stop the entire job before copying; use -f to remove and overwrite those files.
The -p value may include arguments, e.g. -p './postprocess.sh --mode fast'.
Quotes within that value preserve spaces in arguments; shell operators and
variable expansion are not interpreted.
Requires alien.py on PATH and an authenticated GRID session.

Example:
    ./copyFromAlien.py -d folders.txt -l o2_cosmics_full.root,o2_tfidinfo.root \
        -o /existing/output -n 4 -p ./postprocess.sh
"""

import argparse
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path, PurePosixPath
import shutil
import shlex
import subprocess
import sys


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def split_list(value):
    return [item.strip() for item in value.split(",") if item.strip()]


def normalize_folder(value):
    if value.startswith("alien://"):
        value = value[len("alien://"):]
    if ".." in value.split("/"):
        raise ValueError("parent directory components are not allowed: " + value)
    relative = value.lstrip("/")
    if not relative or str(PurePosixPath(relative)) == ".":
        raise ValueError("an empty or root GRID folder is not allowed")
    return "/" + str(PurePosixPath(relative))


def process_folder(folder, files, output, alien, postprocess, overwrite):
    destination = output / folder.lstrip("/")
    # Also reject existing symlinks that would redirect output outside its root.
    destination.resolve().relative_to(output)
    destination.mkdir(parents=True, exist_ok=True)
    for name in files:
        target = destination / name
        if os.path.lexists(target):
            if not overwrite:
                raise FileExistsError("{} already exists; use -f to overwrite existing files".format(target))
            target.unlink()
        print("[{}] Copying {}".format(folder, name), flush=True)
        subprocess.run(
            [alien, "cp", folder + "/" + name, "file:" + str(destination) + "/"],
            check=True,
        )
    if postprocess:
        print("[{}] Running {}".format(folder, shlex.join(postprocess)), flush=True)
        subprocess.run(postprocess, cwd=str(destination), check=True)
    print("[{}] Done".format(folder), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", required=True, metavar="PATHS",
                        help="comma-separated GRID folders or an existing text file (one folder per line)")
    parser.add_argument("-l", "--files", required=True, metavar="FILES",
                        help="comma-separated file names to copy from each folder")
    parser.add_argument("-f", "--force", action="store_true",
                        help="remove existing requested files before copying")
    parser.add_argument("-p", metavar="SCRIPT",
                        help="quoted executable and arguments; run in each copied folder")
    parser.add_argument("-o", default=".", metavar="OUTPUT",
                        help="existing destination directory (default: current directory)")
    parser.add_argument("-n", type=positive_int, default=1, metavar="NPAR",
                        help="maximum folders processed in parallel (default: 1)")
    if len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args()

    output = Path(args.o).expanduser().resolve()
    if not output.is_dir():
        parser.error("output directory must already exist: " + str(output))
    try:
        paths_file = Path(args.d).expanduser()
        if paths_file.is_file():
            paths = [line.strip() for line in paths_file.read_text().splitlines() if line.strip()]
        else:
            paths = split_list(args.d)
        folders = list(dict.fromkeys(normalize_folder(path) for path in paths))
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if not folders:
        parser.error("no GRID folders provided")
    files = list(dict.fromkeys(split_list(args.files)))
    if not files or any(name in (".", "..") or "/" in name for name in files):
        parser.error("--files must contain file names without directory components")
    # Preflight every folder before submitting work, so a conflict anywhere
    # prevents any transfers or postprocessing elsewhere.
    try:
        for folder in folders:
            destination = output / folder.lstrip("/")
            destination.resolve().relative_to(output)
            for name in files:
                target = destination / name
                if os.path.lexists(target):
                    if not args.force:
                        parser.error("{} already exists; use -f to overwrite existing files".format(target))
                    if target.is_dir() and not target.is_symlink():
                        parser.error("requested file is an existing directory: " + str(target))
    except (OSError, ValueError) as error:
        parser.error(str(error))
    alien = shutil.which("alien.py")
    if not alien:
        parser.error("alien.py was not found on PATH")
    postprocess = None
    if args.p is not None:
        try:
            postprocess = shlex.split(args.p)
        except ValueError as error:
            parser.error("invalid postprocessing command: " + str(error))
        if not postprocess:
            parser.error("postprocessing command must not be empty")
        # Resolve before changing directories, including bare names in the cwd.
        candidate = Path(postprocess[0]).expanduser()
        executable = str(candidate.resolve()) if candidate.is_file() else shutil.which(postprocess[0])
        if not executable or not os.access(executable, os.X_OK):
            parser.error("postprocessing script is not executable or was not found: " + postprocess[0])
        postprocess[0] = executable

    failures = 0
    with ThreadPoolExecutor(max_workers=args.n) as pool:
        pending = {pool.submit(process_folder, folder, files, output, alien, postprocess, args.force): folder
                   for folder in folders}
        for future in as_completed(pending):
            try:
                future.result()
            except (OSError, ValueError, subprocess.CalledProcessError) as error:
                failures += 1
                print("[{}] FAILED: {}".format(pending[future], error), file=sys.stderr, flush=True)
    print("Completed: {} succeeded, {} failed".format(len(folders) - failures, failures))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
