import argparse
import csv
import glob
import os
import sys
from typing import List, Tuple

from schrodinger.structure import Structure, StructureReader, _StructureAtom
from schrodinger.structutils.analyze import evaluate_asl


class Args:
    def __init__(self, dir: str, resids: str, output: str, cutoff: float = 5.0) -> None:
        self.dir = dir
        self.cutoff = cutoff
        self.output = output
        tokens = [resid.split(":") for resid in resids.split(",")]
        self.resids = sorted((ch.strip(), int(num)) for ch, num in tokens)


def get_contact_residues(structure: Structure, cutoff: float) -> List[Tuple[str, int]]:
    idxs: List[int] = evaluate_asl(structure, f"within {cutoff} (ligand)")
    atoms: List[_StructureAtom] = [structure.atom[i] for i in idxs]
    return sorted(set((atom.chain, atom.resnum) for atom in atoms))


def parse_args(argv: List[str]) -> Tuple[Args, str]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dir", help="Folder containing Maestro (*.mae[gz]) files")
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=5,
        help="Contact cutoff. Defaults to %(default)s A.",
    )
    parser.add_argument(
        "-r",
        "--residues",
        required=True,
        help="""
        Comma-separated list of residue ids (CHAIN:NUMBER) to report, e.g.,
        'A:81,B:205'.
        """,
    )
    parser.add_argument("-o", "--output", required=True, help="Output CSV filename.")

    opts = parser.parse_args(argv)
    return Args(**vars(opts)), opts.dir


def main(argv: List[str]) -> None:
    opts, dir = parse_args(argv)

    with open(opts.output, "w") as io:
        csvfile = csv.writer(io)

        headers = ["NAME", "INDEX", "IFDSCORE", "DGBIND"]
        headers.extend([f"{ch.upper()}:{resnum}" for ch, resnum in opts.resids])
        csvfile.writerow(headers)

        maefiles = sorted(glob.glob(os.path.join(dir, "*.mae*")))
        for maefile in maefiles:
            basename = os.path.basename(maefile)
            basename = basename.replace(".mae", "").replace("-out", "")
            with StructureReader(maefile) as reader:
                for i, structure in enumerate(reader):
                    resids = get_contact_residues(structure, opts.cutoff)
                    mask = [1 if resid in resids else 0 for resid in opts.resids]

                    name: str = structure.property.get("s_m_title", basename)
                    ifdscore: float = structure.property["r_psp_IFDScore"]
                    dg_bind: float = structure.property["r_psp_MMGBSA_dG_Bind"]
                    csvfile.writerow([name, i, ifdscore, dg_bind] + mask)


if __name__ == "__main__":
    main(sys.argv[1:])