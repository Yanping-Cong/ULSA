import numpy as np
import healpy as hp
import h5py
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import sys
import os
#from ULSA.spectral_index_fitting.spectral_index_analysis_for_constant import reborn
from ULSA.spectral_index_fitting.spectral_index_frequency_dependent import freq_dependent_index
class smooth(object):
    def __init__(self, nside, v, index_type,I_E_form,using_raw_diffuse,using_default_params,input_spectral_index=None, v_file_dir = None, beta_1 = 0.7,v_1 = 1.):
        self.nside = nside
        self.v = v
        self.index_type = index_type
        self.I_E_form = I_E_form
        self.using_raw_diffuse = using_raw_diffuse
        self.input_spectral_index = input_spectral_index
        self.using_default_params = using_default_params
        self.beta_1 = beta_1
        self.v_1 = v_1
        self.v_file_dir = v_file_dir
        _path = os.path.dirname(os.path.abspath(__file__))
        _path = os.path.split(_path)[0]
        _path = os.path.split(_path)[0]
        self.file_dir = _path +'/obs_sky_data'

    def I_E(self, v):
        if self.I_E_form == 'seiffert':
            result = 1.2*(v*1e-3/1.0)**-2.58

        if self.I_E_form == 'seiffert_freq_depend':
            #suppose
            beta_1 = self.beta_1
            v_1 = self.v_1
            beta0 = -2.58
            beta = beta0 + beta_1*np.exp(-v/v_1)
            result = 1.2*(v*1e-3/1.0)**-2.58

        return result

    def free_free_408(self,nside):
        #in galactic coordinate
        with h5py.File(self.file_dir + '/Free_free_emission/T_ff_nu_408MHz.hdf5','r') as f:
            free_free_G = f['free_free'][:]
        free_free_G = hp.ud_grade(free_free_G,nside)
        free_free_G = self.masked_smoothing(free_free_G,5.0)
        return free_free_G

    def freq_dependence_index_minus_I_E(self,freq):
        if self.input_spectral_index != None:
            beta0, beta_1, v_1 = self.input_spectral_index[0],self.input_spectral_index[1],self.input_spectral_index[2]
        if self.input_spectral_index == None:
            if self.using_default_params == True:
                beta0 = -2.51
                beta_1 = self.beta_1;v_1 = self.v_1
                beta = beta0 + beta_1 * np.exp(-freq/v_1)
            if self.using_default_params == False:
                f = freq_dependent_index(freq,self.beta_1,self.v_1,self.v_file_dir)
                beta0 = f.calculate_index(self.nside)
                beta_1 = self.beta_1;v_1 = self.v_1

            beta = beta0 + beta_1 * np.exp(-freq/v_1)
        index_ = beta
        return index_

    def produce_data(self):
        
        data_freq = 408.
        ##the based data of 408MHz coming from HS14 data set
        data_diffuse = hp.read_map(self.file_dir + '/Haslam/haslam408_dsds_Remazeilles2014.fits')
        data_diffuse = hp.ud_grade(data_diffuse, self.nside)
            
        
        if self.index_type == 'constant_index' or 'direction_dependent_index':
            #no need to interpolating to v MHz, because we just do fitting in 408MHz
            diffuse_x = data_diffuse - self.I_E(data_freq)

            mask = np.where(diffuse_x <0)[0]
            diffuse_x[mask] = np.nan

        if self.index_type == 'freq_dependent_index':
            _index = self.freq_dependence_index_minus_I_E(self.v)
            data_diffuse = data_diffuse - self.I_E(data_freq)
                
            diffuse_x = np.multiply(data_diffuse, (self.v/data_freq)**_index)
            mask = np.where(diffuse_x <0)[0]
            diffuse_x[mask] = np.nan

        return diffuse_x

    def nan_helper(self,y):
        """Helper to handle indices and logical indices of NaNs.

        Input:
            - y, 1d numpy array with possible NaNs
        Output:
            - nans, logical indices of NaNs
            - index, a function, with signature indices= index(logical_indices),
              to convert logical indices of NaNs to 'equivalent' indices
        Example:
            >>> # linear interpolation of NaNs
            >>> nans, x= nan_helper(y)
            >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
        """

        return np.isnan(y), lambda z: z.nonzero()[0]

    def masked_smoothing(self, U, rad=5.0):     
        V=U.copy()
        V[U!=U]=0
        VV=hp.smoothing(V, fwhm=np.radians(rad))    
        W=0*U.copy()+1
        W[U!=U]=0
        WW=hp.smoothing(W, fwhm=np.radians(rad))    
        return VV/WW

    def mask_spur_loopI(self, arr):
        nside = self.nside
        l1,b1 = hp.pix2ang(nside,100,lonlat = True)
        l2,b2 = hp.pix2ang(nside,101, lonlat = True)
        #print (l1,l2,b1,b2)
        delt_l = np.abs(l2-l1) 
        #nside = self.nside
        mask = np.zeros_like(arr,dtype = np.bool)
        for l in np.arange(-30,45,0.01):
            for b in np.arange(20,90,delt_l * 0.01):
                pix_number = hp.ang2pix(nside,l,b,lonlat = True)
                mask[pix_number] = True
        arr[mask] = np.nan
        return arr


    def add_5(self):
        diffuse = self.produce_data()

        #remove small scale terbulance, the galactic coordinate list below
        #CAS-a 111.73,-2.13
        #Cygnus 80,4
        #Cen-A 309.51589573409,19.41727341133
        #SMC 302.8084 -44.3277
        #LMC 280.4652 -32.8884
        #Vela 263.9390 -03.3683
        nside = self.nside

        a = hp.get_all_neighbours(nside,111.73,-2.13,lonlat = True)
        b = (hp.get_all_neighbours(nside,80,4,lonlat = True))
        bb = (hp.get_all_neighbours(nside,80,3,lonlat = True))
        c = (hp.get_all_neighbours(nside,309.51589573409,19.41727341133,lonlat = True))
        d = (hp.get_all_neighbours(nside,302.8084,-44.3277,lonlat = True))
        e = (hp.get_all_neighbours(nside,280.4652,-32.8884,lonlat = True))
        ee = (hp.get_all_neighbours(nside,278.4652,-32.8884,lonlat = True))
        f = (hp.get_all_neighbours(nside,263.9390,-3.3683,lonlat = True))
        ff = (hp.get_all_neighbours(nside,263.9390,-2.3683,lonlat = True))
        fff = (hp.get_all_neighbours(nside,263.9390,-1.3683,lonlat = True))
        ffff = (hp.get_all_neighbours(nside,263.9390,0,lonlat = True))
        fffff = (hp.get_all_neighbours(nside,262.9390,-2,lonlat = True))
        ffffff = (hp.get_all_neighbours(nside,264.9390,-2,lonlat = True))
        fffffff = (hp.get_all_neighbours(nside,264.9390,-1,lonlat = True))
        ffffffff = (hp.get_all_neighbours(nside,262.9390,-1,lonlat = True))
        g = (hp.get_all_neighbours(nside,0,0,lonlat = True))
        total = list(a) + list(b) + list(bb) + list(c)+list(d)+list(e)+list(ee)+list(f)+list(ff)+list(fff)+list(ffff)+list(fffff)+list(ffffff)+list(fffffff)+list(ffffffff)+list(g)        

        a = (hp.ang2pix(nside,111.73,-2.13,lonlat = True))
        b = (hp.ang2pix(nside,80,4,lonlat = True))
        c = (hp.ang2pix(nside,309.51589573409,19.41727341133,lonlat = True))
        d = (hp.ang2pix(nside,302.8084,-44.3277,lonlat = True))
        e = (hp.ang2pix(nside,280.4652,-32.8884,lonlat = True))
        f = (hp.ang2pix(nside,263.9390,-3.3683,lonlat = True))
        g = (hp.ang2pix(nside,0,0,lonlat = True))
        total2 = [a,b,c,d,e,f,g]
        pass_pix_num = np.array(total+total2)

 
        # first mask spur and loop I, then smoothing, orthwise, the spur and loop I will reflect the pixel value around spur and loop I 
        diffuse = self.mask_spur_loopI(diffuse) 
        diffuse[pass_pix_num] = np.nan


        #if self.using_raw_diffuse == True:
        #    smooth_diffuse = diffuse
        #    smooth_diffuse = self.mask_spur_loopI(smooth_diffuse)
        #    smooth_diffuse[pass_pix_num] = np.nan
            

        if self.using_raw_diffuse == False or True:
            smooth_diffuse = self.masked_smoothing(diffuse,rad=5.0)
            smooth_diffuse = smooth_diffuse - self.free_free_408(self.nside)
            #mask the value of (raw_map - free - I_E) < 1
            smooth_diffuse[np.where(smooth_diffuse<1.)[0]] = np.nan
            # mask again in case the after smooth will change the mask value 
            smooth_diffuse = self.mask_spur_loopI(smooth_diffuse)
            smooth_diffuse[pass_pix_num] = np.nan

        result = []
        for index in range(smooth_diffuse.size):
            ll,bb = hp.pix2ang(self.nside, index, lonlat = True)
            if smooth_diffuse[index] == np.nan:
                pass
                #print ('pass_pix_number',index)
            if smooth_diffuse[index] != np.nan:
                pix_value = smooth_diffuse[index]
                if ~np.isnan(pix_value):
                    result.append([ll,bb,pix_value])
        result = np.array(result)
        return result
























 
