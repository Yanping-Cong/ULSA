Tutorial
========

.. note::

   This is intended to be a tutorial for the user of *ULSA* package, who will
   just use the already presented tasks in the package to do some data analysis.


.. contents::

konwledge of the input parameters
-------------------------------




This parameter can be used when you are using the fuction in the class::

  * v int: frequency in MHz of the input skymap. 
  * nside int: the NSIDE of the skymap in healpix mode. 
  * clumping_factor float: the clumping factor influnce the absorption, the value set to one in our model. 
  * index_type str: ('constant_index_minus_I_E', 'freq_dependence_index_minus_I_E', 'pixel_dependence_index_minus_I_E'), one of them can be choose as different type of spectral index one need to consider.
  * distance kpc: the integration max distance for galaxy absorption. 
  * test: normally False, if True will do some test programm. 
  * emi_form str: ('exp','sech') the distribution form of emissivity, normally choosen 'exponantial'. 
  * I_E_form str: ('seiffert'), the form of extragalactic component except for CMB. 
  * R0_R1_equal bool: in this paper we fixed R0 equals to R1 in emissivity. 
  * using_raw_diffue bool: the input data for fitting parameter of emissivity, if True the data will be smoothed by Gaussian function. 
  * using_default_params: if True, using the default spectral index value, if False calculate the spectral index value with the code, otherwise, one can simply input the spectral index to variable of using_default_params. 
  * params_408: if the input of params_408 == [0.,0.,0.,0.,0.], the code will fit the parameters of emissivity in 408Mhz, or one can simply input the parameters of some other fitting result to params_408, if you input nothing, the code will take the default parameters.
  * critical_dis: if True, calculate the critial distance (time consuming), otherwise False.
  * output_absorp_free_skymap: if True, output the absorption free sky map in input frequency. 


Run the code
----------------
an example.py example for calculating constant spectral index situation:
    >>> from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
    >>> f = absorption_JRZ(v = v, nside = nside, clumping_factor = 1., index_type = 'constant_index_minus_I_E', distance = dist, test = False, emi_form  = 'exp',I_E_form = 'seiffert',R0_R1_equal=True,using_raw_diffuse = False)
    >>> sky_map_in_healpix = f.mpi()


Single process run
^^^^^^^^^^^^^^^^^^

If you do not have an MPI environment installed, or you just want a single
process run, just do ::

   $ python example.py


If you want to submit and run the pipeline in the background, do like ::

   $ nohup python dir/example/examle.py &> output.txt &

Multiple processes run
^^^^^^^^^^^^^^^^^^^^^^

To run the pipeline in parallel and distributed maner on a cluster using
multiple processes, you can do something like ::

   $ mpiexec -n N python example.py 

or (in case *script.py* isn't in you working directory) ::

   $ mpiexec -n N python dir/example/example.py

If you want to submit and run the pipeline in the background on several nodes,
for example, *node2*, *node3*, *node4*, do like ::

   $ nohup mpiexec -n N -host node2,node3,node4 --map-by node python dir/example/example.py &> output.txt &

.. note::

   In the above commands, **N** is the number of processes you want to run!


products and intermediate results
------------------------------------------

script.py products and intermediate results will be in the running directory in hdf5 file or an array store in your return variable.


