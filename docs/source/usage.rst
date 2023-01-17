Usage
=====

This page shows how to use :py:mod:`polypharm` both programmatically or from
command line.

From Python
-----------

First, import the module:

.. code-block:: python

    import polypharm as ppm

As stated in :doc:`overview`, there are four main steps in the designed
workflow, where each one is executed by a single function:

- Induced-fit docking (IFD) [:py:func:`polypharm.cross_dock`]
- Docking rescore [:py:func:`polypharm.rescore_docking`]
- Report [:py:func:`polypharm.report`]
- Ranking [:py:func:`polypharm.rank_poses` and :py:func:`polypharm.rank_molecules`]

Let's first import utility modules and set a few variables shared across
the steps:

.. code-block:: python

    import glob
    import os

    # required to run docking and binding energy calculation
    os.environ["SCHRODINGER_PATH"] = "/path/to/schrodinger"

    workdir = "/path/to/folder"
    # these are the proteins found in the examples
    raw_resids = {
        "kv1.5_open_cc": "A:436,A:439,A:440,A:443,A:502,A:510,B:436,B:439,B:440",
        "nav1.5_6lqa": "B:1462,B:1466,B:1760,B:1767",
        "task1_6rv3": "A:126,A:171,A:194,A:198,A:199,A:234,A:235,A:236,A:238",
        "kv1.5_open_cc_sp": "A:436,A:439,A:440,A:443,A:502,A:510,B:436,B:439",
    }
    # transform comma-separated strings to a list of strings to avoid typing
    resid_map = {k: v.split(",") for k, v in raw_resids.items()}

The first step is run the induced-fit docking (IFD) for each ligand in
every protein:

.. code-block:: python

    prot_files = glob.glob(os.path.join(workdir, "proteins", "*.mae*"))
    lig_files = glob.glob(os.path.join(workdir, "molecules", "*.mae*"))
    ppm.cross_dock(
        prot_files=prot_files,
        lig_files=lig_files,
        bs_residues=resid_map,
        workdir=os.path.join(workdir, "ifd"),
        glide_cpus=2,
        prime_cpus=2,
        tasks=10,
    )

This will run 10 parallel IFD calculations for every protein/ligand
combination, where Glide and Prime use 2 cores each.

Then, the obtained protein-ligand complexes are rescored using MM/GBSA:

.. code-block:: python

    ifd_files = glob.glob(os.path.join(workdir, "ifd", "**", "*", "*-out.maegz"))
    ppm.rescore_docking(
        ifd_files,
        workdir=os.path.join(workdir, "mmgbsa"),
        cpus=2,
        tasks=10,
    )

Similar to :py:func:`polypharm.cross_dock` function, the above will run
10 parallel MM/GBSA calculations for each output of IFD, where Prime
uses 2 cores.

Once the calculations are done, call the :py:func:`polypharm.report`
function to analyze the output files and generate a
:py:class:`pandas.DataFame` containing the relevant information for
further processing:

.. code-block:: python

    maefiles = glob.glob(os.path.join(workdir, "mmgbsa", "**", "*-out.maegz"))
    results = ppm.report(maefiles, resid_map, contact_cutoff=6, tasks=10)

This will analyze the output files of the MM/GBSA calculations in
parallel using a distance cutoff of 6 Å to detect contacts between the
protein's binding site residues and the ligand atoms.

The fourth step involves sorting of the docking poses and ranking the
molecules, which can be done by one single call:

.. code-block:: python

    criteria = [
        ppm.RankingCriterion.NORMALIZED_CONTACTS,
        ppm.RankingCriterion.TOTAL_SCORE,
    ]
    ranked_results = ppm.rank_molecules(results, criteria)

This sorts the docking poses according to the normalized contacts and
total score (see :py:class:`polypharm.RankingCriterion`) and then ranks
the molecule across the multiple proteins such that those with higher
chance to bind to the targets are ranked first.

The dataframes can be examined further or write to CSV files
(recommended) as follows:

.. code-block:: python

    results.to_csv("/path/to/csv")

It is recommended to check the `pandas's user guide
<https://pandas.pydata.org/docs/user_guide/index.html>`_ to know how to
manipulate (e.g., sorting, filtering, etc.) the dataframes.

Command line
------------

The same workflow can be run via command line as follows.

First, let's try to execute the ``polypharm`` program:

.. code-block:: bash

    $ python -m polypharm
    usage: polypharm [-h] {dock,rescore,report,rank} ...
    polypharm: error: the following arguments are required: {dock,rescore,report,rank}

Use the ``-h`` option to print the help information:

.. code-block:: bash

    $ python -m polypharm
    usage: polypharm [-h] {dock,rescore,report,rank} ...

    Run a stage of the structure-based drug design workflow for polypharmacology.

    positional arguments:
    {dock,rescore,report,rank}
        dock                Run induced-fit cross-docking
        rescore             Run MM/GBSA for the cross-docking output
        report              Generate a report for MM/GBSA output
        rank                Rank molecules by the given criteria across the multiple receptors

    optional arguments:
    -h, --help            show this help message and exit

This is a recurring option in all the tools, so be sure to run every
command with the ``-h`` option at least once to read their
documentation.

It can be seen that there are four commands, each related to the four
main stages shown above, where one can follow the same workflow almost
one to the code. The main difference is that the dataframes are
place by CSV files.

Thus, all the steps can be run sequentially:

.. code-block:: bash

    $ vim resids.txt # write residues into a file
    $ mkdir ifd && cd ifd
    $ python -m polypharm dock -p ../proteins -r ../resids.txt -t 10 ../ligands
    $ cd ..
    $ mkdir mmgbsa && cd mmgbsa
    $ python -m polypharm rescore -t 10 ../ifd
    $ cd ..
    $ python -m polypharm report -o report.csv -c 6 -t 10 mmgbsa
    $ python -m polypharm rank -o rank.csv report.csv

Note that either files or directories can be passed to any command,
where the appropiate files will be searched within the directories
(e.g., ``*-out.maegz`` will be search in the MM/GBSA folder).

Residues can be given as strings or as a plain-text file (``'/path/to/resid'``) containing one entry per line, e.g.,

.. code-block:: text

    kv1.5_open_cc     A:436,A:439,A:440,A:443,A:502,A:510,B:436,B:439,B:440
    nav1.5_6lqa       B:1462,B:1466,B:1760,B:1767
    task1_6rv3        A:126,A:171,A:194,A:198,A:199,A:234,A:235,A:236,A:238
    kv1.5_open_cc_sp  A:436,A:439,A:440,A:443,A:502,A:510,B:436,B:439
