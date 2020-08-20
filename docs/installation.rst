.. _sec-installation:

============
Installation
============

``iBioGen`` requires Python >= 3.5. Installation is facilitated by the conda package
management system.

1. Download `miniconda <https://conda.io/miniconda.html>`_ and run the installer: ``bash Miniconda*``
2. Create a separate `conda environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_ to install iBioGen into:

.. code:: bash

    conda create -n iBioGen
    conda activate iBioGen

3. Install:

.. code:: bash

    conda install -c conda-forge -c iBioGen iBioGen

4. Test:

.. code:: bash

   iBioGen -v

Installation issues can be reported on the `iBioGen github <https://github.com/isaacovercast/iBioGen>`_.
