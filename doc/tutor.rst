Tutorial
========

.. note::

   This is intended to be a tutorial for using of *ULSA* package, who will
   just use the already presented functions in the package to do some simulation.


How to use the code
-------------------

After you sucessfully install the package of ULSA, one can simply do in ipython ::

    >>> from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
    >>> f = absorption_JRZ(v, nside, clumping_factor, index_type, distance,emi_form,I_E_form,R0_R1_equal,using_raw_diffuse,test=False, using_default_params=True, params_408 = np.array([71.19, 4.23, 0.03, 0.47, 0.77]),critical_dis=False,output_absorp_free_skymap=False,beta_1=0.7,v_1 = 1.0)
    >>> sky_map_in_healpix = f.mpi()

This parameter is useful when you are using the fuction:

  * v int: frequency in MHz of the output skymap. 
  * nside int: the NSIDE decide the resolution of the skymap in healpix mode. 
  * clumping_factor float: the clumping factor influnce the absorption, we use the value of 1. in our paper. 
  * index_type str: ('constant_index_minus_I_E', 'freq_dependence_index_minus_I_E', 'pixel_dependence_index_minus_I_E'), one of them can be choose as different model of spectral index.
  * distance kpc: the max integrated distance for galaxy when considering the free-free absorption. 
  * test: normally False, if True one can do some test programm in the position under test==True. 
  * emi_form str: ('exp','sech') the distribution form using in emissivity, 'exponantial' or 'sech^2' will be select, one can change the form of emissivity by need. 
  * I_E_form str:  (’seiffert’,’seiffert freq depend’), the form of extragalactic component minus CMB, if index type == ’constant index minus I E’ or ’pixel dependence index minus I E’, accordingly I E form should　choose ’seiffert’. if index type == ’freq dependence index minus I E’, accordingly I E form should choose ’seiffert freq depend’ . 
  * R0_R1_equal bool:  choosing True in this paper, R0 and R1 is the parameters of emissivity. 
  * using_raw_diffue bool:  the input data for fitting parameter of emissivity, if True the data will be smoothed by Gaussian kernel. 
  * using_default_params bool:  if True, using the default spectral index value, if False calculate the spectral index value with the code, otherwise, one can simply input the spectral index to variable of using default params. 
  * params_408 list:  if the input of params 408 == [0.,0.,0.,0.,0.], the code will fit the parameters of emissivity in　408Mhz, or one can simply input the parameters of some other fitting result to params 408, By default, if you　do nothing, the code will take the default parameters.
  * critical_dis bool: if True, calculating the half brightness distance, otherwise False.
  * output_absorp_free_skymap bool:  if True, output the absorption free sky map at frequency v. 

.. note::

   where, the “v” indicates the frequency of the output sky map. Users can output    sky map at different frequencies according to their needs, or use a loop to output    sky map at a series of frequencies. By default, the output sky map coordinates are    galactic coordinate system, users can convert to another coordinate system by the    “change_coord” function, as described in the following example. 
   The “nside” parameter controls the resolution of the output sky map. In general,    the pixel number of the sky map is equal to 12 * nside * nside. The “clumping    factor” amplify the free-free absorption effect, and in our final model it is    appropriate to determine that the value of the clumping factor is 1. 
   “distance” is set to the galactic maximum integration scale, and in general, 50    kpc is the absorption scale that would be sufficient to cover the entire galaxy. For    normal users, the "test" parameter is set to False, and for developers who want to    modify the code, they can do some custom tests under the code block of test ==    True. “emi_form” determines the distribution of emissivity, so we have two    options, one is an exponential distribution, the other is an exponential distribution    multiplied by the square of a hyperbolic sine function. 
   In the original emissivity distribution, R0 and R1 are two different parameters, and    it is proved that R0 and R1 equaling to each other is an optimal choice. The “using_raw_diffuse” parameter determines whether to smooth the data before fitting the emissivity function, emissivity is a smooth changed distribution with R and Z. The selection of data smoothing can better remove the negative effects of fitting caused by small-scale structure. The remaining parameters are the default parameters, we recommend using the default values, for developers, according to the above parameter description, can flexible to change parameter by they need.

Some examples for different need
-------------------------------------

1. example one: Users want to output a sky map at 1 Mhz with a NSIDE of 64 and a spectral index in the form of
a constant spectral index. they can choose the following parameter setting::

    >>> (v = 1, nside = 64, clumping factor = 1., index type = ’constant index minus I E’, distance = 50,test = False, emi form = ’exp’,I E form = ’seiffert’,R0 R1 equal=True,using raw diffuse = False,using default params = True,critical dis = False,output absorp free skymap = False)

2. example two: Users want to output a sky map at 1 Mhz with a NSIDE of 64 and a spectral index in the form of
a frequency dependent spectral index. they can choose the following parameter setting::

    >>> (v = 1, nside = 64, clumping factor = 1., index type = ’freq dependence index minus I E’, distance = 50, test = False, emi form = ’exp’,I E form = ’seiffert freq depend’,R0 R1 equal = True,using raw diffuse = False,using default params = True,critical dis = False,output absorp free skymap = False)

3. example three: Users want to output a sky map at 1 Mhz with a NSIDE of 64 and a spectral index in the form
of a pixel dependent spectral index. they can choose the following parameter setting::

    >>> (v = 1, nside = 64, clumping factor = 1., index type = ’pixel dependence index minus I E’, distance = 50, test = False, emi form = ’exp’,I E form = ’seiffert’,R0 R1 equal = True,using raw diffuse = False,using default params = True,critical dis = False,output absorp free skymap = False)

Change the coordinate of the output sky_map
--------------------------------------------------

if one want to change the coordinate from 'Galactic' to other coordinate, one can simplely using function "change_coord"::

    >>> from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
    >>> f = absorption_JRZ(v, nside, clumping_factor, index_type, distance,emi_form,I_E_form,R0_R1_equal,using_raw_diffuse,test=False, using_default_params=True, params_408 = np.array([71.19, 4.23, 0.03, 0.47, 0.77]),critical_dis=False,output_absorp_free_skymap=False,beta_1=0.7,v_1 = 1.0)
    >>> sky_map_in_healpix = f.mpi()
    >>> f.change_coord(sky_map_in_healpix,["G","C"])

where ["G","C"], the first character is the coordinate system of sky_map_in_healpix, second character is the coordinate system of the output map. As in HEALPIX, allowed　coordinate systems are 'G' (galactic), 'E' (ecliptic) or 'C' (equatorial)

Run the code
----------------

an example.py example for calculating constant spectral index condition is under ULSA/example/example.py ::

    $ from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
    $ import healpy as hp
    $ import numpy as np
    $ import matplotlib.pyplot as plt
    $ sky_map_list = []
    $ for v in [1,3,10]:
    $     f = absorption_JRZ(v = v, nside = 64, clumping factor = 1., index type = ’pixel dependence index minus I E’, distance = 50, test = False, emi form = ’exp’,I E form = ’seiffert’,R0 R1 equal = True,using raw diffuse = False,using default params = True,critical dis = False,output absorp free skymap = False)
    $     sky_map_list.append(f.mpi())
    $ # we got a list of sky_map with frequency from 1Mhz to 10Mhz with step 0.1Mhz.
    $ # then plot the data using mollview
    $ plt.figure(1)
    $ for sky_map in sky_map_list
    $     hp.mollview(np.log10(sky_map),cmap = plt.cm.jet)
    $     plt.show() # or plt.savefig('xxx.eps',format='eps')

.. note ::

   All the used observation data is from website, they all locate in the dir of /obs_sky_data, if there are new observation data in low frequency, you can change the input data by replacing or adding the data under the dir of /obs_sky_data/*


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


