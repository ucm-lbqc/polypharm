Installation
============

This page describes how to download, install and use the basic
functionality of **polypharm**.

Requirements
------------

This package requires Python 3.7 or greater, so be sure to have a
working installation.

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
<https://github.com/maurobedoya/polypharm>`_.

Uncompress the source code, change directory, and install it by:

.. code-block:: bash

    $ unzip polypharm.zip
    $ cd polypharm
    $ python setup.py install .
