Installation
============

Python version
--------------
Both python2 and python3 can run in a single process，If you want to run in parallel, it depends on whether the installed MPI  package is supported by Python2 or Python3. 

Prerequisites
-------------

For the installation and proper work of ``ULSA``, the following packages are
required:

    * h5py_\ , Pythonic interface to the HDF5 binary data format;
    * numpy_\ , Base N-dimensional array package for Python;
    * scipy_\ , The fundamental package for scientific computing with Python;
    * matplotlib_\ , A python 2D plotting library;
    * caput_\ , Cluster Astronomical Python Utilities;
    * NE2001_\, the free electron distribution in galaxy;
    * mpi4py_\, This package provides Python bindings for the Message Passing Interface (MPI) standard;
.. note:: ``ULSA`` can work without MPI support, in which case, only a single
   process is invoked, but in order to process large amounts of data in parallel
   and distributed manner, mpi4py_ is needed.

Installation guide
------------------
first you need to install the NE2001 under the dir of ULSA/NE2001, you shoud into the dir of /NE2001_4python/src.NE2001 and run make .so, then it will produce an libNE2001.so link, and to test if you are success in install, you should go to the dir of bin_NE2001 and run fortran_python.py, if it is success, mean you do the right thing,then you should remember that  location of libNE2001.s0, and replace the location in ULSA/ULSA/sky_map/produce_absorbed_sky_map.py line 39.

.. note:: _Haslam 408MHz: https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_get.cfm
.. note:: _LWA: http://lda10g.alliance.unm.edu/LWA1LowFrequencySkySurvey/
.. note:: Guzman 45MHz: with hdf5 Form, and in galaxy coordinate

As for respect the copyright of data, the observation data under the dir of ULSA/obs_sky_data/ you should download them by youself, we just giving the data link.

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


finally, for example, if you want to run the code in dir of ULSA/example, you should copy all the file under NE2001_4python/bin_NE2001 to ULSA/example, because it is the input parameter of  NE2001.
.. note:: we update the NE2001 as NE2001_4python and make it faster and can produce a link using in python.

.. _h5py: http://www.h5py.org/
.. _healpy: https://pypi.python.org/pypi/healpy
.. _pyephem: http://rhodesmill.org/pyephem/
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org
.. _matplotlib: http://matplotlib.org
.. _caput: https://github.com/zuoshifan/caput/tree/zuo/develop
.. _mpi4py: http://mpi4py.readthedocs.io/en/stable/
.. _NE2001: http://hosting.astro.cornell.edu/~cordes/NE2001/
