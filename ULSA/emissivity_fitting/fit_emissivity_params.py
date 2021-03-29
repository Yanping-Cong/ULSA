#!/usr/bin/env python
# coding: utf-8
#import matplotlib
#matplotlib.use('Agg')
import numpy as np
from numpy import sin,cos,pi
from scipy.integrate import quad
import matplotlib.pyplot as plt
import scipy.constants as C
import healpy as hp
import h5py
import scipy.optimize as optimize
from scipy.integrate import quad

from pylab import cm
import time

from ULSA.emissivity_fitting.produce_data_for_fitting import smooth
#from Params.I_E_term.I_E_equation import I_E
#from Params.interpolate_sky.interpolate_sky_map import produce_index
class free_free(object):

    def __init__(self, v, nside, index_type, dist, using_raw_diffuse,using_default_params,params_408,input_spectral_index=None,v_file_dir=None,emi_form = 'exp'):
        self.v = v#20Mhz frequency in hz
        self.nside = nside       
        self.index_type = index_type
        self.dist = dist
        if self.index_type == 'constant_index':
            self.I_E_form = 'seiffert'
        if self.index_type == 'direction_dependent_index':
            self.I_E_form = 'seiffert'
        if self.index_type == 'freq_dependent_index':
            self.I_E_form = 'seiffert_freq_depend' 
        self.emi_form =  emi_form
        self.R0_R1_equal = True
        self.input_spectral_index = input_spectral_index
        self.v_file_dir = v_file_dir
        self.using_raw_diffuse = using_raw_diffuse
        self.using_default_params = using_default_params
        #self.params_408 = np.array([71.19, 4.23, 0.03, 0.47, 0.77])
        self.params_408 = params_408
    
    def produce_xyz(self):
        #v in Mhz
        result = smooth(self.nside,self.v, self.index_type,self.I_E_form,self.using_raw_diffuse, self.using_default_params, self.input_spectral_index, self.v_file_dir).add_5()
        return result

    def sech2(self,x):
        return np.square(2/(np.exp(x) + np.exp(-x)))

    def fun_for_curvefit_R0R2(self, xyz, A_v, R_0, alpha, Z_0, gamma, R_2 = 0.1):
        #produce for input for curve_fit
        beta = 1
        R_2 = R_2
        result = []
        for l,b in xyz:
            self.l = l * np.pi / 180.0
            self.b = b * np.pi / 180.0

            def fun_inte(r):
                R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
                Z = r * np.sin(self.b)
                if self.emi_form == 'exp':
                    emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_0)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
                if self.emi_form == 'sech':
                    emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_0)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma)

                #get rid of square 
                return emissivity

            result.append(quad(fun_inte, 0, self.dist)[0])
        return np.array(result)

    def curve_fit(self):
        #if self.R0_R1_equal == False:
        if True:
            #A_v, R_0,R_2=0.1, alpha, Z_0, gamma
            guess = [80,5,2,1.6,1]

            func = self.fun_for_curvefit_R0R2
            xyz = self.produce_xyz()
            #print 'xyz.shape',xyz.shape
            beta_ = -2.51 + 0.7 * np.exp(-self.v/1.0)
            A_upper_limit = 10* 15 * (self.v/408.)**beta_
            #print 'A_upper_limit',A_upper_limit
            params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess, bounds=(np.array([0,1e-5,1e-5,1e-5,1e-5]),np.array([A_upper_limit,5,3.1,2,3.1])), method='trf')
        #print ('params_408',params)
        with h5py.File(str(self.v)+'Mhz_fitted_param.hdf5','w') as f:
            f.create_dataset('params',data = params)
            f.create_dataset('v',data = self.v)
            f.create_dataset('pcov', data = pcov)
        return params
 
    def params(self):
        if self.index_type == 'constant_index' or 'direction_dependent_index' or 'freq_dependent_index':

            if self.params_408.all() == 0 or 0.:
                try:
                    with h5py.File('./'+str(self.v)+'Mhz_fitted_param.hdf5','r') as f:
                        abcz0 = f['params'][:]
                except:
                    if int(self.v) == int(408):
                        print ('fitting the params of synchrotron emissivity in 408MHz') 
                        abcz0 = self.curve_fit()
                    else:
                        raise 'fitting params must at frequency of 408MHz for constant ,pixel dependent and freq dependent spectral index situation'
            else:
                
                abcz0 = self.params_408

        return abcz0



