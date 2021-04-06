from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
import numpy as np
#input_spectral_index = np.array([2.5])
import time 
print ('time_begin',time.ctime())
#v =1.
nside = 2**6
dist = 50.
if __name__ == '__main__':
    #cla = absorption_JRZ(v = v, nside = nside,index_type = 'constant_index_minus_I_E', distance = dist,using_raw_diffuse = False,using_default_params=True,input_spectral_index=input_spectral_index)
    for v in np.arange(0.1,1,0.1):
    #for v in [1.,3.,10.]:
        v = round(v,1)
        cla = absorption_JRZ(v = v, nside = nside,index_type = 'constant_index', distance = dist,using_raw_diffuse = False,using_default_params=True,output_absorp_free_skymap=True,critical_dis = False)
        cla.mpi()
        print ('time_end',time.ctime())
