ULSA.emissivity_fitting.produce_data_for_fitting.smooth
=======================================================
class::

   class ULSA.emissivity_fitting.produce_data_for_fitting.smooth(object):

method:

.. py:function:: __init__(self, nside, v, index_type, I_E_form, using_raw_diffuse,using_default_params,beta_1=0.7,v_1=1.)
   
   initial parameter function

   :param int nside: The Nside value one choose in healpix mode, must be 2^N.
   :param v float: The frequency of output sky map
   :param str index_type: ('constant_index_minus_I_E', 'freq_dependence_index_minus_I_E', 'pixel_dependence_index_minus_I_E'), one of them can be choose as different type of spectral index one need to consider.
   :param str I_E_form: ('seiffert'), the form of extragalactic component except for CMB.
   :param bool using_raw_diffuse: if True, using the raw input data for fitting process without smoothing.
   :param bool using_default_params: if True, using the default spectral index value, if False calculate the spectral index value with the code, otherwise, one can simply input the spectral index to variable of using_default_params.
   :param float beta_1: for frequency dependence index condition, beta = beta0 + beta_1(v/v_1).
   :param float v_1: for frequency dependence index condition, beta = beta0 + beta_1(v/v_1).
.. py:function:: add_5()

   return the (l, b, pixel_value) array for fitting the params of emissivity.

   :return: (longtitude, latitude, pixel_value) in healpix mode with Nside was set
   :rtype: np.array

   
   
   


   
   
