from genericpath import isdir
import glob
import os
import re
import sys
import argparse
from typing import Any, Callable, Dict, List, Tuple

import polypharm


def command_report(args: Dict[str, Any]) -> None:
    csvfile: str = args["csvfile"]
    del args["csvfile"]
    df = polypharm.report(**args)
    df.to_csv(csvfile, index=False)


def command_ifd(args: Dict[str, Any]) -> None:
    args["prot_files"] = expand_files(args["prot_files"])
    args["lig_files"] = expand_files(args["lig_files"])
    polypharm.run_ifd_cross(**args)


def expand_files(paths: List[str]) -> List[str]:
    files: List[str] = []
    for path in paths:
        if os.path.isdir(path):
            files.extend(glob.glob(os.path.join(path, "*.mae*")))
        else:
            files.append(path)
    return files


def parse_residues(raw: str) -> Dict[str, List[str]]:
    if ":" not in raw and os.path.exists(raw):
        with open(raw) as io:
            raw = ";".join(
                [re.sub(r"\s+", "=", line.strip()) for line in io.readlines()]
            )
    pairs = [token.split("=") for token in raw.split(";")]
    pairs = [(prot, raw_resids.split(",")) for prot, raw_resids in pairs]
    return dict(pairs)


def parse_args(
    argv: List[str],
) -> Tuple[Callable[[Dict[str, Any]], None], Dict[str, Any]]:
    main_parser = argparse.ArgumentParser(
        description=__doc__,
        prog="polypharm",
    )
    subparsers = main_parser.add_subparsers(required=True)

    parser = subparsers.add_parser(
        "report", help="Generate a report for MM/GBSA output"
    )
    parser.add_argument(
        "output_dir",
        metavar="DIR",
        help="""
        Directory containing MM/GBSA output files (*-out.maegz) from a cross-docking 
        experiment.
        """,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="csvfile",
        metavar="FILE",
        required=True,
        help="CSV file to write the report to.",
    )
    parser.add_argument(
        "-r",
        "--residues",
        dest="bs_residues",
        metavar="LIST",
        type=parse_residues,
        required=True,
        help="""
        Binding site residues to report. Must be specified as a 
        semicolon-separated list of <protein name>=<comma-separated list of 
        residues>, e.g., 'prot1=A:12,A:64,B:3;prot2=A:34,A:82,C:23,D:24'.
        Alternatively, a file can be given containing one entry per line (the
        '=' symbol may be replaced by a space or tab).
        """,
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        dest="contact_cutoff",
        metavar="FLOAT",
        type=float,
        default=5,
        help="""
        Residue contact cutoff in Å. Selected residues are reported to be within the 
        given cutoff. Defaults to %(default)s Å.
        """,
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_false",
        dest="use_existing",
        help="""
        Residue contact cutoff in Å. Selected residues are reported to be within the 
        given cutoff. Defaults to %(default)s Å.
        """,
    )
    parser.add_argument(
        "-t",
        "--tasks",
        type=int,
        default=1,
        help="Number of parallel tasks to run. Defaults to %(default)s.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Do not print progress to standard output.",
    )
    parser.set_defaults(cmd=command_report)

    parser = subparsers.add_parser("ifd", help="Run induced-fit cross-docking")
    parser.add_argument(
        "lig_files",
        metavar="FILE_OR_DIR",
        nargs=argparse.ONE_OR_MORE,
        help="""
        Maestro files (*.mae[gz]) containing prepared ligands. Directories can 
        also be passed, where *.mae[gz] files will be listed.
        """,
    )
    parser.add_argument(
        "-p",
        "--protein",
        dest="prot_files",
        metavar="FILE_OR_DIR",
        action="append",
        help="""
        Maestro files (*.mae[gz]) containing prepared proteins. Directories can 
        also be passed, where *.mae[gz] files will be listed.
        """,
    )
    parser.add_argument(
        "-r",
        "--residues",
        dest="bs_residues",
        metavar="LIST",
        type=parse_residues,
        required=True,
        help="""
        Binding site residues to report. Must be specified as a 
        semicolon-separated list of <protein name>=<comma-separated list of 
        residues>, e.g., 'prot1=A:12,A:64,B:3;prot2=A:34,A:82,C:23,D:24'.
        Alternatively, a file can be given containing one entry per line (the
        '=' symbol may be replaced by a space or tab).
        """,
    )
    parser.add_argument(
        "-t",
        "--tasks",
        type=int,
        default=1,
        help="Number of parallel tasks to run. Defaults to %(default)s.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Do not print progress to standard output.",
    )
    parser.add_argument(
        "--glide-cpus",
        type=int,
        default=1,
        help="Number of processors to be used by Glide.",
    )
    parser.add_argument(
        "--prime-cpus",
        type=int,
        default=1,
        help="Number of processors to be used by Prime MM/GBSA.",
    )
    parser.set_defaults(cmd=command_ifd)

    args = main_parser.parse_args(argv)
    cmd: Callable[[Dict[str, Any]], None] = args.cmd
    del args.cmd

    return cmd, vars(args)


def main(argv: List[str]) -> None:
    func, args = parse_args(argv)
    func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
