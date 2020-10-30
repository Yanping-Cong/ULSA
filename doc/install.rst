Installation
============

Python version
--------------
Both python2 and python3 can run in a single process，If you want to run in parallel, it depends on whether the installed MPI  package is supported by Python2 or Python3. 

Prerequisites
-------------

For the installation and proper work of ``LFSM``, the following packages are
required:

    * h5py_\ , Pythonic interface to the HDF5 binary data format;
    * numpy_\ , Base N-dimensional array package for Python;
    * scipy_\ , The fundamental package for scientific computing with Python;
    * matplotlib_\ , A python 2D plotting library;
    * caput_\ , Cluster Astronomical Python Utilities;
    * NE2001_\, the free electron distribution in galaxy;
    * mpi4py_\, This package provides Python bindings for the Message Passing Interface (MPI) standard;
.. note:: ``LFSM`` can work without MPI support, in which case, only a single
   process is invoked, but in order to process large amounts of data in parallel
   and distributed manner, mpi4py_ is needed.

Installation guide
------------------

After you have successfully installed the prerequisites, do the following.

First clone this package ::

    $ git clone https://github.com/Yanping-Cong/ULSA

Then change to the top directory of this package, install it by the usual
methods, either the standard ::

    $ python setup.py install [--user]

or to develop the package ::

    $ python setup.py develop [--user]

It should also be installable directly with `pip` using the command ::

    $ pip install [-e] git+https://github.com/Yanping-Cong/ULSA.git




.. _h5py: http://www.h5py.org/
.. _healpy: https://pypi.python.org/pypi/healpy
.. _pyephem: http://rhodesmill.org/pyephem/
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org
.. _matplotlib: http://matplotlib.org
.. _caput: https://github.com/zuoshifan/caput/tree/zuo/develop
.. _mpi4py: http://mpi4py.readthedocs.io/en/stable/
.. _NE2001: http://hosting.astro.cornell.edu/~cordes/NE2001/
