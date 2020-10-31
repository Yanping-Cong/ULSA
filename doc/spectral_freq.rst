ULSA.spectral_index_fitting.spectral_index_analysis_for_direction_dependence
============================================================================ 
class::

   class ULSA.spectral_index_fitting.spectral_index_frequency_dependent.freq_dependent_index(object):

method:

.. py:function:: __init__(self,freq,beta_1,v_1)

   :param float freq: spectral index value in this frequency will be output.
   :param float beta_1: beta = beta0 + beta_1*exp(freq/v_1)
   :param float v_1: beta = beta0 + beta_1*exp(freq/v_1)
  
   :return: none

.. py:function:: combined_index(downgrade_to)

   return the spectral index sky map for spatial-variance spectral index condition.

   :param int downgrade_to: The Nside value one choose in healpix mode, must be 2^N.
   :return: the spectral index value
   :rtype: float
   :raises TypeError: if the downgrade_to is not 2^N.

   
   
   


   
   
