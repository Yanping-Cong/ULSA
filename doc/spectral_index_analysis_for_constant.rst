ULSA.spectral_index_fitting.spectral_index_analysis_for_constant
================================================================

class::

   class ULSA.spectral_index_fitting.spectral_index_constant.constant_index(object):

method:

.. py:function:: __init__(self)
   
   initial parameter function
   :return: none
.. py:function:: calculate_index(downgrade_to)

   return the spectral index value for constant spectral index condition.

   :param int downgrade_to: The Nside value one choose in healpix mode, must be 2^N.
   :return: the spectral index value
   :rtype: float
   :raises TypeError: if the downgrade_to is not 2^N

   
   
   


   
   
