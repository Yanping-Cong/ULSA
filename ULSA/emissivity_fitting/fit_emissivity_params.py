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

    def __init__(self, v, nside, index_type, dist, using_raw_diffuse,using_default_params,input_spectral_index=None,v_file_dir=None,params_408=np.array([71.19, 4.23, 0.03, 0.47, 0.77])
):
        self.v = v#20Mhz frequency in hz
        self.nside = nside       
        self.index_type = index_type
        self.dist = dist
        if self.index_type == 'constant_index_minus_I_E':
            self.I_E_form = 'seiffert'
        if self.index_type == 'pixel_dependence_index_minus_I_E':
            self.I_E_form = 'seiffert'
        if self.index_type == 'freq_dependence_index_minus_I_E':
            self.I_E_form = 'seiffert_freq_depend' 
        self.emi_form = 'exp'
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
                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_0)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
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
            #params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess,sigma = xyz[:,2], bounds=(np.array([0,1e-5,-3.1,1e-5,-3.1]),np.array([1e10,100,3.1,20,3.1])), method='trf')
            #params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess, bounds=(np.array([0,1e-5,-3.1,1e-5,-3.1]),np.array([1e10,100,3.1,20,3.1])), method='trf')
            beta_ = -2.49 + 0.7 * np.exp(-self.v/1.0)
            A_upper_limit = 10* 15 * (self.v/408.)**beta_
            print 'A_upper_limit',A_upper_limit
            params, pcov = optimize.curve_fit(func, xyz[:,:2], xyz[:,2], guess, bounds=(np.array([0,1e-5,1e-5,1e-5,1e-5]),np.array([A_upper_limit,5,3.1,2,3.1])), method='trf')

        with h5py.File(str(self.v)+'Mhz_fitted_param.hdf5','w') as f:
            f.create_dataset('params',data = params)
            f.create_dataset('v',data = self.v)
            f.create_dataset('pcov', data = pcov)
        return params
 
    def params(self):
        if self.index_type == 'constant_index_minus_I_E' or 'pixel_dependence_index_minus_I_E':

            if self.params_408.all() == 0 or 0.:
                try:
                    with h5py.File('./'+str(self.v)+'Mhz_fitted_param.hdf5','r') as f:
                        abcz0 = f['params'][:]
                except:
                    if int(self.v) == int(408):
                        abcz0 = self.curve_fit()
                    else:
                        raise 'fitting params must at frequency of 408MHz for constant and pixel dependence spectral index situation'
            else:
                abcz0 = self.params_408

        if self.index_type == 'freq_dependence_index_minus_I_E':
            abcz0 = self.curve_fit()
        return abcz0



