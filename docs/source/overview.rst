Overview
========

:py:mod:`polypharm` is an open source library written in Python to
perform IFD and MM/GBSA calculations on different targets using a
polypharmacological approach.

Features
--------

It allows performing Induced Fit Docking (IFD) and binding free energy
calculations using the MM/GBSA method available in the Schrödinger suite
by means of a wrapper library. The protocol performs the calculations on
different targets/proteins and considering the experimentally relevant
residues creates a scoring function to evaluate the best candidates with
activity on the different targets, i.e. it proposes the molecules with
the best polypharmacological activity.


Workflow
--------

The workflow can be divided in four main stages:

- :ref:`docking`
- :ref:`rescore`
- :ref:`report`
- :ref:`ranking`

.. _docking:

Cross docking
~~~~~~~~~~~~~

.. warning::

    The initial structures must be previously prepared, with added
    hydrogens, correct bond orders and charges, and assignment of 
    protonation states to the desired pH. Verify that everything is 
    correct in the initial .mae and .maegz files.

An example of the initial setup using a Jupyter notebook can be 
found here:
`usage.ipynb <https://github.com/ucm-lbqc/polypharm/blob/main/examples/usage.ipynb>`_
    
First, import the libraries, including the Polypharm library, and 
set the environment variable for the Schrodinger suite. The names 
of the input and output folders have to be specified. The "LIGAND_DIR" 
folder must contain all the ligands to be tested (in mae or maegz format) 
and the "PROTEIN_DIR" folder all the proteins to be studied 
(in mae or maegz format). If more than one binding site has to be 
specified for the same protein, the file must be duplicated and renamed 
to identify the results. This name must be specified in the form:
    
.. code-block:: python

    raw_residues = {"proteinA_siteA": "A:436,A:439,A:440",
                    "proteinA_siteB": "E:400,E:401,E:402"}
Note that it is possible to use any name with alphanumeric characters, but we 
suggest using simple and short names for the ligands, e.g. 1a, 1b, 1c, ... etc.    

Then the number of core processes to be used in each calculation and the number 
of concurrent processes to be run must be specified. In the following code:

.. code-block:: python

    ppm.cross_dock(
        prot_files=prot_files,
        lig_files=lig_files,
        bs_residues=BS_RESIDS,
        workdir=IFD_DIR,
        glide_cpus=2,
        prime_cpus=2,
        tasks=10,
    )

The Glide and Prime stages use two processors and run 10 calculations
simultaneously, so 20 processor threads are required. In this procedure, 
cross-docking is performed, i.e. docking with the IFD protocol of each ligand 
in each protein and defined binding site and a maximum of 10 poses are generated 
for each ligand in each defined site.

.. _rescore:

Docking rescoring
~~~~~~~~~~~~~~~~~
    
In this step, a rescoring of the docking poses is performed using the MMGBSA 
method to obtain an approximation of the binding free energy for all the generated poses.
Rescoring is defined with a simple function, as follows:

.. code-block:: python

    ppm.rescore_docking(
        glob.glob(os.path.join(IFD_DIR, "**", "*", "*-out.maegz")),
        workdir=MMGBSA_DIR,
        cpus=2,
        tasks=10,
    )

.. _report:

Docking reporting
~~~~~~~~~~~~~~~~~

The report creates a data frame containing all the information generated in the IFD and 
MMGBSA steps. The data frame also contains the binding site residues with which each 
output pose interacts.
The report can be safely run as follows:

.. code-block:: python

    import subprocess
    try:
        maefiles = glob.glob(os.path.join(MMGBSA_DIR, "**", "*-out.maegz"))
        results = ppm.report(maefiles, BS_RESIDS, CONTACT_CUTOFF, tasks=10)
    except subprocess.CalledProcessError as ex:
        print(" ".join(ex.cmd))
        print(ex.stdout.decode())
    results.head(10)

.. _ranking:

Ranking
~~~~~~~

To rank the poses and visualize the energies for each protein/binding site and 
also write the data to a .csv file, the following code can be used:

.. code-block:: python

    ranked_results = ppm.rank_poses(results, RANK_CRITERIA)
    for protein, pdf in ranked_results.groupby("PROTEIN"):
        pdf = pdf.loc[pdf.groupby("NAME")["RANK"].idxmin()]
        pdf = pdf.sort_values("RANK")
        pdf = pdf.copy()
        pdf["RANK"] = list(range(1,len(pdf) + 1))
        pdf = pdf.dropna(axis=1)[["NAME", "INDEX", "DGBIND", "INT", "DGBIND_NORM", "INT_NORM", "NORMT", "RANK"]]
        pdf.to_csv(f"{protein}_rank_poses.csv")
        print(protein)
        display(pdf)
Finally, the molecules are ranked to find those with the greatest polypharmacological 
potential, i.e. those with the lowest global ranking ("GLOBAL RANK") among all the proteins 
and binding sites studied. The following code can be used to obtain the ranking of the 
molecules:

.. code-block:: python

    cross_results = ppm.rank_molecules(results, RANK_CRITERIA)
    cross_results.to_csv("ranking.csv", index=True)
    cross_results
