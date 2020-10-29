from ULSA.sky_map.produce_healpix_sky_map import absorption_JRZ
v =1.
nside = 2**4
dist = 50.
if __name__ == '__main__':
    cla = absorption_JRZ(v = v, nside = nside, clumping_factor = 1., index_type = 'constant_index_minus_I_E', distance = dist, test = False, emi_form  = 'exp',I_E_form = 'seiffert',R0_R1_equal=True,using_raw_diffuse = False)
    cla.mpi()
