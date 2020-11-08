from ULSA.sky_map.produce_absorbed_sky_map import absorption_JRZ
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

def produce_sky_map():
    sky_map_list = []
    # initial params setting
    nside = 2**6
    dist = 50.
    for v in range(1,10,0.1):
        f = absorption_JRZ(v = v, nside = nside, clumping_factor = 1., index_type = 'constant_index_minus_I_E', distance = dist, test = False, emi_form  = 'exp',I_E_form = 'seiffert',R0_R1_equal=True,using_raw_diffuse = False,using_default_params=False,critical_dis = False,output_absorp_free_skymap = False)
        sky_map_list.append(f.mpi())
    # we got a list of sky_map with frequency from 1Mhz to 10Mhz with step 0.1Mhz.
    # then plot the data using mollview
    return sky_map_list
 
def plot():
    sky_map_list = produce_sky_map()
    plt.figure(1)
    for sky_map in sky_map_list
        hp.mollview(np.log10(sky_map),cmap = plt.cm.jet)
        plt.show() # or plt.savefig('xxx.eps',format='eps')

plot()
