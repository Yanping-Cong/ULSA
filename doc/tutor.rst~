Tutorial
========

.. note::

   This is intended to be a tutorial for using of *ULSA* package, who will
   just use the already presented functions in the package to do some simulation.


How to use the code
-------------------

After you sucessfully install the package of ULSA, one can simply do in ipython ::

    >>> from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
    >>> f = absorption_JRZ(v, nside, index_type, distance, using_raw_diffuse, using_default_params,critical_dis,output_absorp_free_skymap)
    >>> sky_map_in_healpix = f.mpi()

This parameter is useful when you are using the fuction:
  	* v (float): frequency in MHz for the output map;
	* nside (int): the Healpix NSIDE value for the sky map;
	* index_type (str): ('constant\_index', 'freq\_dependent\_index', 'direction\_dependent\_index'), spectral index modeling option;
	* distance kpc: maximum integration distance along line of sight, default is set to 50 kpc;
	* using_raw_diffue (bool): If False, the data will be smoothed by Gaussian kernel, otherwise use the raw data;
	* v_file_dir (dict):  a dictionary structure used to specify additional input map data. Specify the frequency of the map by the dictionary key, and the relative path of the map data as dictionary value. The  input sky map file should be in HDF5 format, such as \{XX:"/dir/xxx.hdf5"\}. If None, the spectral index is calculated with the existing data.
	* using_default_params (bool): if the input is a bool type, use the default spectral index value in the code if it is set to True, or re-calculate the spectral index value if it is set to False.  
	* input_spectral_index (array): one can specify the spectral index value by putting in an array containing the spectral index map in the direction dependent case or containing one element for constant or frequency dependent cases.}
	* params_408 (list): The emissivity model parameters ($[A,R_0,\alpha,Z_1,\gamma]$) obtained by fitting at 408 MHz.  If this parameter is omitted, the values given in Table\ref{table_params} will be used as defaults. One can also specify these parameters directly by putting in the values, or force the code to re-fit by setting it to [0.,0.,0.,0.,0.] 
	* critical_di (bool): if True, calculating the half brightness distance, otherwise this is not calculated. 
	* output_absorp_free_skymap (bool): if True, output an addional absorption-free sky map at frequency {\bf v}.  


.. note::

   where, the “v” indicates the frequency of the output sky map. Users can output sky map at different frequencies according to their needs, or use a loop to output sky map at a series of frequencies. By default, the output sky map coordinates are galactic coordinate system, users can convert to another coordinate system by the “change_coord” function, as described in the following example. The “nside” parameter controls the resolution of the output sky map. In general, the pixel number of the sky map is equal to 12 * nside * nside.  “distance” is set to the galactic maximum integration scale, and in general, 50 kpc is the absorption scale that would be sufficient to cover the entire galaxy. The “using_raw_diffuse” parameter determines whether to smooth the data before fitting the emissivity function, emissivity is a smooth changed distribution with R and Z. The selection of data smoothing can better remove the negative effects of fitting caused by small-scale structure. 

Some examples for different need
-------------------------------------

1. example one: Users want to output a sky map at 1 Mhz with a NSIDE of 64 and a spectral index in the form of
a constant spectral index. they can choose the following parameter setting::

    >>> (v = 1, nside = 64, index_type = ’constant_index_minus’, distance = 50, using_raw_diffuse = False,using_default_params = True,critical_dis = False,output_absorp_free_skymap = False)

2. example two: Users want to output a sky map at 1 Mhz with a NSIDE of 64 and a spectral index in the form of
a frequency dependent spectral index. they can choose the following parameter setting::

    >>> (v = 1, nside = 64, index_type = ’freq_dependent_index’, distance = 50, using_raw_diffuse = False,using_default_params = True,critical_dis = False,output_absorp_free_skymap = False)

3. example three: Users want to output a sky map at 1 Mhz with a NSIDE of 64 and a spectral index in the form
of a direction dependent spectral index. they can choose the following parameter setting::

    >>> (v = 1, nside = 64, index_type = ’direction_dependent_index’, distance = 50, using_raw_diffuse = False,using_default_params = True,critical_dis = False,output_absorp_free_skymap = False)

Change the coordinate of the output sky_map
--------------------------------------------------

if one want to change the coordinate from 'Galactic' to other coordinate, one can simplely using function "change_coord"::

    >>> from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
    >>> f = absorption_JRZ(v, nside, index_type, distance,using_raw_diffuse, using_default_params=True,input_spectral_index = None, params_408 = np.array([71.19, 4.23, 0.03, 0.47, 0.77]),critical_dis=False,output_absorp_free_skymap=False)
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
    $ #calculating the sky map from 1Mhz to 10Mhz
    $ for v in [1,10,1]:
    $     #nside in healpix mode
    $     nside = 2**6
    $     f = absorption_JRZ(v = v, nside = nside, index_type = 'constant_index', distance = dist,using_raw_diffuse = False,using_default_params=False,input_spectral_index = None, critical_dis = False,output_absorp_free_skymap = False)
    $     sky_map_list.append(f.mpi())
    $ # we got a list of sky_map with frequency from 1Mhz to 10Mhz with step 1Mhz.
    $ # then plot the data using mollview
    $ plt.figure(1)
    $ for sky_map in sky_map_list
    $     hp.mollview(np.log10(sky_map),cmap = plt.cm.jet)
    $     plt.show() # or plt.savefig('xxx.eps',format='eps')

.. note ::

   All the used observation data is from website, they all locate in the dir of /obs_sky_data, if there are new observation data in low frequency, you can change the input data by replacing or adding the data under the dir of /obs_sky_data/*
   By default, the value of the v\_file\_dir = None parameter set to None, and if you want to add new data, it should be passed to the parameter of v\_file\_dir by a dictionary whose key is the frequency of the input sky map and whose value is the path relative /obs\_sky\_data/* .
   The default file format is HDF5, and the default ‘key’ of data in HealPix mode is ‘data’, and the default coordinate should be 'C' (equatorial).\\
   For example, a file of 22MHz\_sky\_map.hdf5 putting under /obs\_sky\_data/22MHz/22MHz\_sky\_map.hdf5, you should give a dict \{22:"/22MHz/22MHz\_sky\_map.hdf5"\} to v\_file\_dir as v\_file\_dir = \{22:"/22MHz/22MHz\_sky\_map.hdf5"\}, the code will adding the new data when calculating the spectral index automatically.


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


