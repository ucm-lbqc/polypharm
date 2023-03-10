Installation
============

This page describes how to download, install and use the basic
functionality of :py:mod:`polypharm`.

Requirements
------------

- `Python <https://python.org>`_ 3.9 or later
- `pandas <https://pandas.pydata.org>`_ 1.4.3 or later
- `Jinja2 <https://palletsprojects.com/p/jinja>`_ 3.1.2 or later

.. warning::

    The main functionality (i.e., IFD and MM/GBSA) does require a
    working `Schrödinger Suite <https://schrodinger.com>`_ installation
    (2018-4 or greater) including the Glide and Prime modules.

Installing using pip
--------------------

`pip <https://pip.pypa.io/en/stable/>`_ is a command-line tool, which
often comes pre-packaged with Python distributions, that can download
packages from the `Python Package Index (PyPI) <https://pypi.org>`_. To
see if it's installed, on Linux or macOS try:

.. code-block:: bash

    $ which pip

Then, the package can be installed globally via:

.. code-block:: bash

    $ python -m pip install polypharm

or locally for the current user:

.. code-block:: bash

    $ python -m pip install --user polypharm

To ensure the package is available, try:

.. code-block:: bash

    $ python -m polypharm -h
    # or
    $ python -c 'import polypharm; print(polypharm.__version__)'

Manual installation
-------------------

The source code including unreleased changes can be also directly
download from the repository at `GitHub
<https://github.com/ucm-lbqc/polypharm>`_.

Uncompress the source code, change directory, and install it by:

.. code-block:: bash

    $ unzip polypharm.zip
    $ cd polypharm
    $ python setup.py install .
