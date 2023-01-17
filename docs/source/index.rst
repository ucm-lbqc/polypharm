.. polypharm documentation master file, created by
   sphinx-quickstart on Mon Jan 16 13:41:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to polypharm's documentation!
=====================================

:py:mod:`polypharm` is an open source library written in Python to
perform IFD and MM/GBSA calculations on different targets using a
polypharmacological approach.

About
-----

:py:mod:`polypharm` was designed and implemented as part of an
investigation regarding blockers of ion channels. Such research aimed at
identifying possible candidates via computational methods for lead
optimization as part of a drug design workflow. The main goal was to
find one molecule that could potentially modulate several ion channels
related to known diseases, approach often referred to as
polypharmacology.

This investigation was developed at the Centro de Bioinformática y
Simulación Molecular (CBSM_), Universidad de Talca, Chile, led by `Prof.
Wendy Gonzalez <wgonzalez_>`_, in collaboration with the Laboratorio de
Bioinformática y Química Computacional (LBQC_), Universidad Católica del
Maule, Chile, with the help of `Mauricio Bedoya <mbedoya_>`_ (PhD) and
`Francisco Adasme <fadasme_>`_ (PhD).

Please refer to the :doc:`overview` for further details about the
features and methodology.

How to use polypharm
--------------------

The package can be installed via the `pip
<https://pip.pypa.io/en/stable/>`_ package manager or download directly
from the repository at GitHub_ by following the instructions at the
:doc:`installation` page.

:py:mod:`polypharm` can be used either programmatically via the
available Python package or from the command line using the distributed
executable of the same name. Please refer to the :doc:`usage` page for
details.

The documentation of the application programming interface (API), that
is, description of functions and classes, can be found at the :doc:`api`
page.

If you find a bug, need a missing feature, or have any question, please
submit your inquiry via the issue tracker at the GitHub_ repository.

Citing
------

If you use :py:mod:`polypharm` in your scientific work, please support
our work by adding a reference to the following article:

   To be added

Contributors
------------

- `Mauricio Bedoya <https://github.com/maurobedoya>`__ - creator,
  maintainer
- `Francisco Adasme <https://github.com/franciscoadasme>`_ - maintainer

License
-------

Licensed under the MIT license, see the separate LICENSE file.

.. _cbsm: https://cbsm.utalca.cl/
.. _fadasme: https://www.researchgate.net/profile/Francisco-Adasme
.. _github: https://github.com/maurobedoya/polypharm
.. _lbqc: https://lbqc.ucm.cl
.. _mbedoya: https://www.researchgate.net/profile/Mauricio-Bedoya-2
.. _wgonzalez: https://www.researchgate.net/profile/Wendy-Gonzalez-9

.. toctree::
   :maxdepth: 2
   :caption: Contents
   :hidden:

   overview
   installation
   usage
   api