ULSA.emissivity_fitting.fit_emissivity_params.free_free
=======================================================
class::

   class ULSA.emissivity_fitting.fit_emissivity_params.free_free(object):

method:

.. py:function:: __init__(self,v,nside,index_type,dist,emi_form,I_E_form,R0_R1_equal,using_raw_diffuse,using_default_params,params_408=np.array([71.19, 4.23, 0.03, 0.47, 0.77]))
   
   initial parameter function

   :param v float: The frequency of output sky map
   :param int nside: The Nside value one choose in healpix mode, must be 2^N.
   :param str index_type: ('constant_index_minus_I_E', 'freq_dependence_index_minus_I_E', 'pixel_dependence_index_minus_I_E'), one of them can be choose as different type of spectral index one need to consider.
   :param int dist: the maximux integrated distance of galaxy, normally setting to 50kpc.
   :param emi_form:  ['exp','sech'] the distribution form of emissivity, normally choosen 'exponantial'.
   :param str I_E_form: ('seiffert'), the form of extragalactic component except for CMB.
   :param bool R0_R1_equal: fixed True
   :param bool using_raw_diffuse: if True, using the raw input data without smoothing.
   :param bool using_default_params: if True, using the default spectral index value, if False calculate the spectral index value with the code, otherwise, one can simply input the spectral index to variable of using_default_params. 
   :param array params_408: if the input of params_408 == [0.,0.,0.,0.,0.], the code will fit the parameters of emissivity in 408Mhz, or one can simply input the parameters of some other fitting result to params_408, if you input nothing, the code will take the default parameters. 
.. py:function:: params()

   return the fitting params of emissivity in array form.

   :return: params in array form $array([A, R_0, \alpha, Z_1, \gamma])$.
   :rtype: np.array

   
   
   


   
   
