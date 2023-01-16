Overview
========

:py:mod:`polypharm` is an open source library written in Python to perform IFD
and MM/GBSA calculations on different targets using a polypharmacological
approach.

Features
--------

It allows performing Induced Fit Docking (IFD) and binding free energy
calculations using the MM/GBSA method available in the Schrödinger suite by means
of a wrapper library. The protocol performs the calculations on different
targets/proteins and considering the experimentally relevant residues creates a
scoring function to evaluate the best candidates with activity on the different
targets, i.e. it proposes the molecules with the best polypharmacological
activity.

**NOTE**: The initial structures must be previously prepared, with added
hydrogens, correct bond orders and charges, and assignment of protonation states
to the desired pH. Verify that everything is correct in the initial .mae and
.maegz files.
