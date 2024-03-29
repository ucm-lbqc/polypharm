import argparse
import glob
import itertools
import os
import re
import sys
from functools import partial
from typing import Any, Callable, Dict, List, Tuple

import pandas as pd
import polypharm


def command_dock(args: Dict[str, Any]) -> None:
    args["workdir"] = os.getcwd()
    polypharm.cross_dock(**args)


def command_rank(args: Dict[str, Any]) -> None:
    df = pd.read_csv(args["csvfile"])
    df = polypharm.rank_molecules(df, args["criteria"])
    df.to_csv(args["output"], index=False)


def command_report(args: Dict[str, Any]) -> None:
    csvfile: str = args["csvfile"]
    del args["csvfile"]
    df = polypharm.report(**args)
    df.to_csv(csvfile, index=False)


def command_rescore(args: Dict[str, Any]) -> None:
    args["workdir"] = os.getcwd()
    polypharm.rescore_docking(**args)


def expand_path(path: str, ext: List[str]) -> List[str]:
    if not os.path.exists(path):
        raise ValueError(f"{path} not found")
    if os.path.isdir(path):
        return list(
            itertools.chain.from_iterable(
                glob.glob(os.path.join(path, f"*{e}")) for e in ext
            )
        )
    elif os.path.splitext(path)[1] not in ext:
        raise ValueError(f"invalid file type: {path}")
    else:
        return [path]


def parse_args(
    argv: List[str],
) -> Tuple[Callable[[Dict[str, Any]], None], Dict[str, Any]]:
    main_parser = argparse.ArgumentParser(
        description="""
        Run a stage of the structure-based drug design workflow for polypharmacology.
        """,
        prog="polypharm",
    )
    subparsers = main_parser.add_subparsers(required=True)

    # Parser for IFD
    parser = subparsers.add_parser("dock", help="Run induced-fit cross-docking")
    parser.add_argument(
        "lig_files",
        metavar="LIGFILE_OR_DIR",
        nargs=argparse.ONE_OR_MORE,
        type=partial(expand_path, ext=[".mae", ".maegz", ".mae.gz"]),
        help="""
        Maestro files (*.mae[gz]) containing prepared ligands. Directories can 
        also be passed, where *.mae[gz] files will be listed.
        """,
    )
    parser.add_argument(
        "-p",
        "--protein",
        dest="prot_files",
        metavar="PROTFILE_OR_DIR",
        type=partial(expand_path, ext=[".mae", ".maegz", ".mae.gz"]),
        action="append",
        required=True,
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
        '=' symbol may be replaced by spaces or tabs).
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
        help="Number of processors to be used by Glide. Defaults to %(default)s.",
    )
    parser.add_argument(
        "--prime-cpus",
        type=int,
        default=1,
        help="Number of processors to be used by Prime. Defaults to %(default)s.",
    )
    parser.set_defaults(cmd=command_dock)

    # Parser for MM/GBSA
    parser = subparsers.add_parser(
        "rescore", help="Run MM/GBSA for the cross-docking output"
    )
    parser.add_argument(
        "ifd_files",
        metavar="IFDFILE_OR_DIR",
        nargs=argparse.ONE_OR_MORE,
        type=partial(expand_path, ext=["-out.mae", "-out.maegz", "-out.mae.gz"]),
        help="""
        Maestro files (*-out.mae[gz]) generated by induced-fit docking (IFD) runs.
        Directories can also be passed, where *-out.mae[gz] files will be listed.
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
        "--prime-cpus",
        dest="cpus",
        type=int,
        default=1,
        help="Number of processors to be used by Prime MM/GBSA.",
    )
    parser.set_defaults(cmd=command_rescore)

    # Parser for reporting
    parser = subparsers.add_parser(
        "report", help="Generate a report for MM/GBSA output"
    )
    parser.add_argument(
        "maefiles",
        metavar="MAEFILE_OR_DIR",
        nargs=argparse.ONE_OR_MORE,
        type=partial(expand_path, ext=["-out.mae", "-out.maegz", "-out.mae.gz"]),
        help="""
        Maestro files (*-out.mae[gz]) generated by MM/GBSA runs. Directories can also
        be passed, where *-out.mae[gz] files will be listed.
        """,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="csvfile",
        metavar="CSVFILE",
        required=True,
        help="CSV file (*.csv) to write the report to.",
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

    # Parser for ranking
    parser = subparsers.add_parser(
        "rank",
        help="Rank molecules by the given criteria across the multiple receptors",
    )
    parser.add_argument(
        "csvfile",
        metavar="CSVFILE",
        help="""
        Comma-separated value file (*.csv) containing the output of the report command.
        """,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        metavar="CSVFILE",
        required=True,
        help="CSV file (*.csv) to write the ranked molecules to.",
    )
    parser.add_argument(
        "-s",
        "--sort-by",
        dest="criteria",
        metavar="SORT_LIST",
        type=parse_ranking_criteria,
        default="normalized_contacts,total_score",
        help="""
        Sort docking poses by the given criteria given as a comma-separated list.
        Valid values are: contacts, normalized_contacts, binding_energy, 
        normalized_binding_energy, total_score. Values are case-insensitive. Defaults to
        '%(default)s'.
        """,
    )
    parser.set_defaults(cmd=command_rank)

    args = main_parser.parse_args(argv)
    cmd: Callable[[Dict[str, Any]], None] = args.cmd
    del args.cmd

    # calling expand_path per argument produces to nested arrays, so flat them
    for name, value in vars(args).items():
        if name.endswith("files"):
            setattr(args, name, list(itertools.chain.from_iterable(value)))

    return cmd, vars(args)


def parse_ranking_criteria(raw: str) -> List[polypharm.RankingCriterion]:
    criteria: List[polypharm.RankingCriterion] = []
    for token in raw.split(","):
        value = getattr(polypharm.RankingCriterion, token.strip().upper(), None)
        if not value:
            raise ValueError(f"invalid sort criterion: {token}")
        criteria.append(value)
    return criteria


def parse_residues(raw: str) -> Dict[str, List[str]]:
    if ":" not in raw and os.path.exists(raw):
        with open(raw) as io:
            raw = ";".join(
                [re.sub(r"\s+", "=", line.strip()) for line in io.readlines()]
            )
    pairs = [token.split("=") for token in raw.split(";")]
    pairs = [(prot, raw_resids.split(",")) for prot, raw_resids in pairs]
    return dict(pairs)


def main(argv: List[str]) -> None:
    func, args = parse_args(argv)
    func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
