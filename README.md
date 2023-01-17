# Polypharm

A Python-based library to perform induced fit docking (IFD) and MM/GBSA
calculations on different targets using a polypharmacological approach.

Refer to the [official documentation](http://polypharm.rtfd.io/) for
details about installation, usage, methodology, and developer interface.

## Installation

The version at the Python Package Index (PyPI) is always the latest
stable release that is relatively bug-free and can be installed via pip:

```shell
pip install polypharm
```

The minimum Python version is 3.9, and requires the pandas and Jinja2
packages. Refer to the
[Installation](http://polypharm.rtfd.io/en/latest/installation.html)
page of the documentation for more details.

**NOTE**: The main functionality (i.e., IFD and MM/GBSA) does require a
working [Schrödinger Suite](https://schrodinger.com) installation (2018-4 or
greater) including the Glide and Prime modules.

## Usage

`polypharm` can be used either programmatically or from the command
line. There is Jupyter Notebook at the examples folder that shows the
common usage. In any case, below is a brief example.

```python
import glob
import os

import polypharm as ppm

# required to run docking and binding energy calculation
os.environ["SCHRODINGER_PATH"] = "/path/to/schrodinger"

# gather input files and configuration
prot_files = glob.glob(os.path.join("proteins", "*.mae*"))
lig_files = glob.glob(os.path.join("molecules", "*.mae*"))
resid_map = {
    "6lqa": ["B:1462", "B:1466", "B:1760", "B:1767"],
    "6rv3": [
        "A:126", "A:171", "A:194", "A:198", "A:199", "A:234", "A:235",
        "A:236", "A:238"
    ],
}
parallel = 10

# 1. Run induced-fit cross docking
ppm.cross_dock(prot_files, lig_files, resid_map, tasks=parallel)
# 2. Rescore generated IFD poses using MM/GBSA
ppm.rescore_docking(
    glob.glob(os.path.join("ifd", "**", "*", "*-out.maegz")),
    tasks=parallel,
)
# 3. Generate a report from MM/GBSA output
maefiles = glob.glob(os.path.join("mmgbsa", "**", "*-out.maegz"))
results = ppm.report(maefiles, resid_map, tasks=parallel)
# 4. Rank molecules by their docking performance
ranked_results = ppm.rank_molecules(results)
```

The same workflow can be performed via command line:

```shell
$ vim resids.txt # write residues into a file
$ mkdir ifd && cd ifd
$ python -m polypharm dock -p ../proteins -r ../resids.txt -t 5 ../ligands
$ cd ..
$ mkdir mmgbsa && cd mmgbsa
$ python -m polypharm rescore -t 5 ../ifd
$ cd ..
$ python -m polypharm report -o report.csv -c 6 -t 5 mmgbsa
$ python -m polypharm rank -o rank.csv report.csv
```

Please refer to the [official documentation](http://polypharm.rtfd.io)
for more information.

## Citing

If you use `polypharm` in your research, please consider citing the
following article:

    To be added

## Contributors

- [Mauricio Bedoya](https://github.com/maurobedoya) - creator,
  maintainer
- [Francisco Adasme](https://github.com/franciscoadasme) - maintainer

## License

Licensed under the MIT license, see the separate LICENSE file.
