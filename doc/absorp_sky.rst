ULSA.sky_map.produce_absorbed_sky_map.absorption_JRZ
===========================================================
class::

   class ULSA.sky_map.produce_absorbed_sky_map.absorption_JRZ(object):

method:

.. py:function:: __init__(self,v,nside,clumping_factor,index_type,distance,emi_form,I_E_form,R0_R1_equal,using_raw_diffuse,test,using_default_params=True,params_408=np.array([71.19,4.23,0.03,0.47,0.77]),critical_dis=False,output_absorp_free_skymap=False,beta_1=0.7,v_1=1.0)
   
   initial parameter function

   :param v float: The frequency of output sky map
   :param int nside: The Nside value one choose in healpix mode, must be 2^N.
   :param int clumping_factor: the clumping factor influnce the absorption, the value set to one in our model. 
   :param str index_type: ('constant_index_minus_I_E', 'freq_dependence_index_minus_I_E', 'pixel_dependence_index_minus_I_E'), one of them can be choose as different type of spectral index one need to consider.
   :param int distance: the maximux integrated distance of galaxy, normally setting to 50kpc.
   :param emi_form:  ['exp','sech'] the distribution form of emissivity, normally choosen 'exponantial'.
   :param str I_E_form: ('seiffert'), the form of extragalactic component except for CMB.
   :param bool R0_R1_equal: fixed True
   :param bool using_raw_diffuse: if True, using the raw input data without smoothing.
   :param test bool test: if True, one can do some test with different parameter using in the code. normally fixed False.
   :param bool using_default_params: if True, using the default spectral index value, if False calculate the spectral index value with the code, otherwise, one can simply input the spectral index to variable of using_default_params. 
   :param array params_408: if the input of params_408 == [0.,0.,0.,0.,0.], the code will fit the parameters of emissivity in 408Mhz, or one can simply input the parameters of some other fitting result to params_408, if you input nothing, the code will take the default parameters.
   :param bool critical_dis: if True, calculate the critial distance (time consuming), otherwise False.
   :param bool output_absorp_free_skymap: if True, output the absorption free sky map in input frequency.
   :param float beta_1: beta = beta0 + beta_1*exp(freq/v_1)
   :param float v_1: beta = beta0 + beta_1*exp(freq/v_1)

.. py:function:: mpi()

   return array([pixel_number, pixel_value]) in healpix mode, shape as (12*nside**2,2). 

   :return: a numpy array shape as [12*Nside^2, 2] will be return one can simply plot in healpy.mollview.
   :rtype: np.array

   
   
   


   
   
