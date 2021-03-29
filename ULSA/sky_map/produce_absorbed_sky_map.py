#!/usr/bin/env python
# coding: utf-8
print ('Version-0.9,separate free-free and synchrotron emission.')
import matplotlib
matplotlib.use('Agg')
import scipy
import h5py
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import astropy.units as u
import healpy as hp

import numpy as np
from numpy import sin,cos,pi
from scipy.integrate import quad
import matplotlib.pyplot as plt
import scipy.constants as C
import healpy as hp
import h5py
import scipy.optimize as optimize
from scipy.integrate import quad

#from matplotlib import cm
from pylab import cm
import time
from caput import mpiutil
from ULSA.emissivity_fitting.fit_emissivity_params import free_free
from ULSA.spectral_index_fitting.spectral_index_constant import constant_index
from ULSA.spectral_index_fitting.spectral_index_frequency_dependent import freq_dependent_index
from ULSA.spectral_index_fitting.spectral_index_direction_dependent import direction_dependent_index
from ULSA.emissivity_fitting.produce_data_for_fitting import smooth
import ctypes as ct
import numpy as np
import os
import sys

# import the dll
libNE2001 = ct.CDLL('/public/home/wufq/congyanping/Software/NE2001_4python/src.NE2001/libNE2001.so')
# max integrated distance (kpc)
dist = 50.


rank = mpiutil.rank
size = mpiutil.size
class absorption_JRZ(object):
    
    def __init__(self, v, nside, index_type,using_raw_diffuse, distance = 50., v_file_dir=None, emi_form='exp',R0_R1_equal=True, using_default_params=True, input_spectral_index = None, params_408 = np.array([43.09932087,3.40820378,0.46230938,1.11894379,1.2344227]),critical_dis=False,output_absorp_free_skymap=False,beta_1=0.7,v_1 = 1.0):
        self.v = v
        self.nside = nside
        self.index_type = index_type
        self.dist = distance
        self.emi_form = emi_form
        self.R0_R1_equal = R0_R1_equal
        self.using_raw_diffuse = using_raw_diffuse
        self.input_spectral_index = input_spectral_index
        self.using_default_params = using_default_params
        self.v_file_dir = v_file_dir
        self.params_408 = params_408
        self.critical_dis = critical_dis
        self.output_absorp_free_skymap = output_absorp_free_skymap
        Te = 8000 
        _path = os.path.dirname(os.path.abspath(__file__))
        _path = os.path.split(_path)[0]
        _path = os.path.split(_path)[0]
        self.file_dir = _path +'/obs_sky_data'
        self.beta_1 = beta_1
        self.v_1 = v_1
        self.file_dir2 = _path +'/obs_sky_data'
        self.file_dir1 = _path +'/ULSA/spectral_index_fitting'
        if self.index_type == 'constant_index':
            self.I_E_form = 'seiffert'
            if rank == 0:

                self.Beta_G_constant = self.constant_index_minus_I_E()
            else:
                self.Beta_G_constant = None
            self.Beta_G_constant = mpiutil.bcast(self.Beta_G_constant, root = 0)
        if self.index_type == 'direction_dependent_index':
            self.I_E_form = 'seiffert'
            if rank == 0:
                self.Beta_G = self.pixel_dependence_index_minus_I_E()
            else:
                self.Beta_G = None
            self.Beta_G = mpiutil.bcast(self.Beta_G, root = 0)
            
        if self.index_type == 'freq_dependent_index':
            self.I_E_form = 'seiffert_freq_depend'
            if rank == 0:
                self.Beta_G_freq = self.freq_dependence_index_minus_I_E(self.v)
            else:
                self.Beta_G_freq = None
            self.Beta_G_freq = mpiutil.bcast(self.Beta_G_freq, root = 0)



    def constant_index_minus_I_E(self):
        if self.input_spectral_index != None:
            beta = self.input_spectral_index[0]
        if self.input_spectral_index == None:
            if self.using_default_params == True:
                beta = -2.51
            if self.using_default_params == False:
                f = constant_index(self.v_file_dir)
                beta = f.calculate_index(self.nside)

        #beta = float(-2.4894)
        index_ = np.ones(12*self.nside**2) * beta
        return index_

    def freq_dependence_index_minus_I_E(self,freq):
        if self.input_spectral_index != None:
            beta0, beta_1, v_1 = self.input_spectral_index[0],self.input_spectral_index[1],self.input_spectral_index[2]
            beta = beta0 + beta_1 * np.exp(-freq/v_1)
        if self.input_spectral_index == None:
            if self.using_default_params == True:
                beta0 = -2.51
                beta_1 = self.beta_1;v_1 = self.v_1
                beta = beta0 + beta_1 * np.exp(-freq/v_1)
            if self.using_default_params == False:
                f = freq_dependent_index(freq, self.beta_1, self.v_1,self.v_file_dir)
                beta0 = f.calculate_index(self.nside)
                beta_1 = self.beta_1;v_1 = self.v_1
                beta = beta0 + beta_1 * np.exp(-freq/v_1)
        index_ = np.ones(12*self.nside**2) * beta
        #print ('index',index_,index_.shape)
        return index_

    def pixel_dependence_index_minus_I_E(self):
        if self.input_spectral_index != None:
            index = self.input_spectral_index
        if self.input_spectral_index == None:
            if self.using_default_params == True:
                try:
                    with h5py.File(self.file_dir1+'/spectral_index_map.hdf5','r') as f:
                        index = f['spectral_index'][:]
                except:
                    f = direction_dependent_index(self.v_file_dir)
                    index = f.combined_index(256)

            if self.using_default_params == False:
                f = direction_dependent_index(self.v_file_dir)
                index = f.combined_index(256)
                plt.figure(1)
                hp.mollview(index,cmap = plt.cm.jet)
                plt.savefig('index_nside_256.png',ppi=600) 
        return hp.ud_grade(index,self.nside)
    
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

    def masked_smoothing(self, U, rad=56/60.):
        V=U.copy()
        V[U!=U]=0
        VV=hp.smoothing(V, fwhm=np.radians(rad))
        W=0*U.copy()+1
        W[U!=U]=0
        WW=hp.smoothing(W, fwhm=np.radians(rad))
        return VV/WW
 
    def free_free_408(self):
        #in galactic coordinate
        with h5py.File(self.file_dir + '/Free_free_emission/T_ff_nu_408MHz.hdf5','r') as f:
            free_free_G = f['free_free'][:]
        free_free_G = self.masked_smoothing(free_free_G,56/60.)
        free_free_G = hp.ud_grade(free_free_G,self.nside)
        return free_free_G

    def diffuse_x(self, freq):

        if self.index_type == 'constant_index':
            #index = self.constant_index_minus_I_E()
            index = self.Beta_G_constant 
        if self.index_type == 'freq_dependent_index':
            index = self.freq_dependence_index_minus_I_E(freq)

        if self.index_type == 'direction_dependent_index':
            #index = self.pixel_dependence_index_minus_I_E()
            index = self.Beta_G 

        
        data_freq = 408.
        #the based data of 408MHz coming from SSM model output
        #data_diffuse, source, catalog = self.SSM_(data_freq)

        ##the based data of 408MHz coming from HS14 data set
        data_diffuse = hp.read_map(self.file_dir2 + '/Haslam/haslam408_dsds_Remazeilles2014.fits')
        #free free having higher resolution, so free_free_408 function have smoothed the skymap to 56 arcmin
        #first ud_grade otherwise too many minus value
        data_diffuse = hp.ud_grade(data_diffuse, self.nside)

        Mask = np.where(data_diffuse - self.free_free_408() - self.I_E(data_freq) <0)[0]
        data_diffuse[Mask] = np.nan
        nans, x= self.nan_helper(data_diffuse)
        data_diffuse[nans]= np.interp(x(nans), x(~nans), data_diffuse[~nans])


        data_diffuse = data_diffuse - self.free_free_408() - self.I_E(data_freq)
        #smooth 408MHz map, 408 is 56 arcmin is the lowest resolution no need to smooth
        #data_diffuse = self.masked_smoothing(data_diffuse)
 
        Mask = np.where(data_diffuse <0)[0]
        data_diffuse[Mask] = np.nan
        nans, x= self.nan_helper(data_diffuse)
        data_diffuse[nans]= np.interp(x(nans), x(~nans), data_diffuse[~nans])
        

        diffuse_x = np.multiply(data_diffuse, (freq/data_freq)**index)

        value = self.I_E(freq)
        if self.output_absorp_free_skymap == True: 
            with h5py.File(str(freq)+'MHz_absorp_free_skymap.hdf5', 'w') as f:
                f.create_dataset('freq', data = freq)
                f.create_dataset('nside',data = self.nside)
                f.create_dataset('diffuse_x', data = diffuse_x)
                f.create_dataset('diffuse_408',data = data_diffuse)
                f.create_dataset('index',data = index)
                f.create_dataset('I_E',data = np.array([value,np.nan]))
        return diffuse_x

    def model_m2(self,l,b,abcz0,R_2=0.1):

        self.l = l * np.pi / 180.0
        self.b = b * np.pi /180.0
        A_v, R_0,alpha, Z_0, gamma = abcz0
        R_1 = R_0.copy()
        R_2 = R_2
        beta = 1
        Beta_G = self.Beta_G_freq
        def fun_inte(r):

            #integrate along the sight direction
            R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
            Z = r * np.sin(self.b)
            emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**Beta_G[0]

            #get rid of square 
            return emissivity 

        return quad(fun_inte, 0, self.dist)[0]

    def model_m4(self,l,b,params_408,pix_number,R_2 = 0.1):
        self.l = l * np.pi / 180.0
        self.b = b * np.pi /180.0
        A_v, R_0,alpha, Z_0, gamma = params_408
        R_1 = R_0.copy()
        R_2 = R_2
        beta = 1
        #Beta_G = f.pixel_dependence_index_minus_I_E()
        #Beta_G = self.pixel_dependence_index_minus_I_E()
        Beta_G = self.Beta_G 
        def fun_inte(r):

            #integrate along the sight direction
            R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
            Z = r * np.sin(self.b)
            emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**Beta_G[pix_number]
            #get rid of square 
            return emissivity 

        return quad(fun_inte, 0, self.dist)[0]

    def model_m5(self,l,b,params_408,R_2 = 0.1):
        self.l = l * np.pi / 180.0
        self.b = b * np.pi /180.0
        A_v, R_0,alpha, Z_0, gamma = params_408
        R_1 = R_0.copy()
        R_2 = R_2
        beta = 1
        #Beta_G = self.constant_index_minus_I_E()
        Beta_G = self.Beta_G_constant 
        def fun_inte(r):

            #integrate along the sight direction
            R = np.sqrt(8.5**2 + (r*np.cos(self.b))**2 -2*8.5*(r*np.cos(self.b))*np.cos(self.l))
            Z = r * np.sin(self.b)
            emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**Beta_G[0]
            #get rid of square 
            return emissivity 
        return quad(fun_inte, 0, self.dist)[0]
    
    def Delt_T(self):
        g = free_free(v = self.v, nside = self.nside,index_type = self.index_type,dist = self.dist,using_raw_diffuse = self.using_raw_diffuse,using_default_params = self.using_default_params,input_spectral_index = self.input_spectral_index, v_file_dir = self.v_file_dir, emi_form = self.emi_form,params_408=self.params_408)
        params = g.params()
        abcz0 = params.copy()
        nside = self.nside
        m = np.zeros(hp.nside2npix(nside))

        I_E = self.I_E(self.v)
        for pix_number in range(0,hp.nside2npix(nside)):
            l,b = hp.pixelfunc.pix2ang(nside, pix_number, nest = False, lonlat = True)
            #the emissivity after absorption and plus the term of extragalaxy
            a = time.time()
            if self.index_type == 'direction_dependent_index':
                pix_value = self.model_m4(l,b,abcz0,pix_number) + I_E
            if self.index_type == 'constant_index':
                pix_value = self.model_m5(l,b,abcz0) + I_E
                    
            if self.index_type == 'freq_dependent_index': 
                pix_value = self.model_m2(l,b,abcz0) + I_E
            m[pix_number] = pix_value
            

            b = time.time()
        #m contains I_E you can see line 279
        with h5py.File(str(self.v)+'m.hdf5','w') as f:
            f.create_dataset('m',data = m)
            f.close()
            
        diffuse_raw = self.diffuse_x(self.v) 
        #m contains I_E you can see line 279
        delt_m = (diffuse_raw + I_E) - m
        return params,delt_m,diffuse_raw,m

    def Fortran2Py_optical_deepth(self, l, b, Te = 8000):
        v = self.v * 1e6 #v in MHz
        v = float(v)
        rad=57.2957795
        #radian per degree
        #distance equals 50kpc
        #dist=50.0
        
        #if self.test == True:
        if False:
            
            step = 0.1
        else:
            step = 0.01
        N =np.int(dist/step)

        nd = ct.pointer( ct.c_int(N) )          # setup the pointer

        em1D = np.arange(0, N, dtype=np.float32)  # setup the N-long

        l_rad = l / rad #now its radian unit
        b_rad = b / rad
        _ = libNE2001.em_los_(nd, ct.pointer( ct.c_float(l_rad) ), ct.pointer( ct.c_float(b_rad) ), ct.pointer( ct.c_float(dist) ), np.ctypeslib.as_ctypes(em1D))
        #EM = pyne2001.get_dm_full(l, b, r)['EM']	
        Tao_mw = 3.28*1e-7 * (Te/1e4)**-1.35 * (v * 1e-9)**-2.1 * em1D
        #print 'Tao_mw',Tao_mw 
        return Tao_mw
 
         
    def integrate_by_hand(self, f, a, b, args = [], dx=0.01):
        #if self.test == True:
        if False:
            dx = 0.1
            step = dx
        else:
            dx = 0.01
            step = dx

        tao = self.Fortran2Py_optical_deepth(args[0], args[1])
        
        i = a 
        s = 0
        I_E = self.I_E(self.v)

        while i <= b:
            #index_ = np.int(i / step - 1)
            index_ = np.int(i / step) - 1
            s += (f(i,args[0],args[1],args[2],args[3]) * np.exp(-tao[index_])) * dx
            i += dx
        #here find the bug
        s = s + I_E*np.exp(-tao[-1])
        return s

    def Quad(self, f, a, b, args = [], dx=0.01):
        #the different to integrate_by_hand is not including I_E
        #if self.test == True:
        if False:
            dx = 0.1
            step = dx
        else:
            dx = 0.01
            step = dx

        tao = self.Fortran2Py_optical_deepth(args[0], args[1])
        
        i = a 
        s = 0
        #take the left value dont take the right value
        while i < b:
            index_ = np.int(i / step - 1)
            s += (f(i,args[0],args[1],args[2],args[3]) * np.exp(-tao[index_])) * dx
            i += dx
        #here find the bug
        s = s 
        return s


    def split_array(self, container, count):
        #return [container[_i::count] for _i in range(count)]
        return np.split(container, count)

    def gaussian(self, x, mu = 8.5, sigma = 1.33333):
        f =  1./np.sqrt(2*np.pi*sigma**2)* np.exp(-(x-mu)**2 / (2*sigma**2))
        return f

    def sech2(self,x):
        return np.square(2/(np.exp(x) + np.exp(-x)))

    def I_E(self, v):
        if self.I_E_form == 'seiffert':
            result = 1.2*(v*1e-3/1.0)**-2.58

        if self.I_E_form == 'seiffert_freq_depend':
            #suppose
            beta_1 = self.beta_1
            v_1 = self.v_1
            beta0 = -2.58
            beta = beta0 + beta_1*np.exp(-v/v_1)
            result = 1.2*(v*1e-3/1.0)**beta

        return result

    def _new(self, r, l, b, delt_m, params):
        if self.R0_R1_equal == True:

            param = params
            A_v = param[0]
            R_0 = param[1]
            R_2 = 0.1
            alpha = param[2]
            R_1 = param[1]
            #beta = param[3]
            beta = 1
            Z_0 = param[3]
            gamma = param[4]

        if self.R0_R1_equal == False:
            param = params
            A_v = param[0]
            R_0 = param[1]
            alpha = param[2]
            R_1 = param[3]
            beta = param[4]
            Z_0 = param[5]
            gamma = param[6]



        r0 = 8.5 

        l_rad = l * np.pi/180.
        b_rad = b * np.pi/180.
        
        R = np.sqrt(8.5**2 + (r*np.cos(b_rad))**2 -2*8.5*(r*np.cos(b_rad))*np.cos(l_rad))
        Z = r * np.sin(b_rad)

        ########ne = (R/(R_0+0.1))**alpha * a * np.exp(-np.abs(Z) * 2/(B+0.1) - 2*(r_1/(20*c + 0.1))**2) + D
        #emissivity = A_v * (R/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
        if self.index_type == 'direction_dependent_index':
            pix_number = hp.ang2pix(self.nside, l, b, lonlat = True)
            if self.emi_form == 'exp':

                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**self.Beta_G[pix_number]
            if self.emi_form == 'sech2':
                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**self.Beta_G[pix_number]

        if self.index_type == 'constant_index':
            if int(self.v) == int(408):
                if self.emi_form == 'exp':
                    emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)
                if self.emi_form == 'sech2':
                    emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma)
            else:
                if self.emi_form == 'exp':

                    emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**self.Beta_G_constant[0]
                if self.emi_form == 'sech2':
                    emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**self.Beta_G_constant[0]
             
        if self.index_type == 'freq_dependent_index':
            if self.emi_form == 'exp':

                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * np.exp(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**self.Beta_G_freq[0]
            if self.emi_form == 'sech2':
                    
                emissivity = A_v * ((R+R_2)/R_0)**alpha * np.exp(-(R/R_1)**beta) * self.sech2(-(np.abs(Z)/Z_0)**gamma)*(self.v/408.)**self.Beta_G_freq[0]
        j_RZ = emissivity #+ delt_m/dist) #* np.exp(-tao[index])
        return j_RZ

    def critical_distance(self,l,b,delt_m,params):
        import scipy.optimize as so
        #bug report : the lower limit is from 0.01 not 0
        value = 0.5 * self.Quad(self._new, 0.01, 50,args=(l,b,delt_m,params)) 
        
        add_value = 0
        Y = []
        for i in np.arange(0.01,6.6,0.01):
            add_value += self.Quad(self._new,i,i+0.01,args=(l,b,delt_m,params))
            Y.append(add_value)
            if add_value - value <= 0:
                sol = i
                continue
        #Y = list(np.abs(Y))
        #container = np.arange(0.01,6.6,0.01)
        #index = Y.index(min(Y))
        #sol = container[index]  
        return sol



    def mpi(self):
        rank = mpiutil.rank
        size = mpiutil.size

        if rank == 0:
            print ('in rank0') 
            params,delt_m,diffuse_raw_,unabsorb  = self.Delt_T()
        else:
            delt_m = None
            params = None
            diffuse_raw_ = None
            unabsorb = None
        #local_delt_m = mpiutil.mpilist(delt_m, method = 'con',comm = MPI.COMM_WORLD)
        local_range = mpiutil.mpirange(0,hp.nside2npix(self.nside))

        delt_m = mpiutil.bcast(delt_m, root = 0)
        params = mpiutil.bcast(params, root = 0)
        diffuse_raw_ = mpiutil.bcast(diffuse_raw_, root = 0)
        unabsorb = mpiutil.bcast(unabsorb, root = 0)
        result_absorb = []
        for pix_number in local_range:
            a = time.time()
            l, b = hp.pix2ang(self.nside, pix_number, nest = False, lonlat = True)
            #if self.test == True:
            if False:
                pix_value =self.integrate_by_hand(self._new, 0.1, dist, args=(l, b, delt_m[pix_number], params)) 
            else:

                pix_value =self.integrate_by_hand(self._new, 0.01, dist, args=(l, b, delt_m[pix_number], params))
                if self.critical_dis == True: 
                    distance = self.critical_distance(l,b,delt_m[pix_number],params)
                    #print 'distance',distance
                else:
                    distance = 0.
                l, b = hp.pix2ang(self.nside, pix_number, nest = False, lonlat = True)
            b = time.time()
            
            #if self.test == True:
            if False:

                result_absorb.append([pix_number, pix_value])
            else:
                result_absorb.append([pix_number, pix_value, distance])
            print ('rest_number',12*self.nside**2 - pix_number)
        #if self.test == True:
        if False:

            result_absorb = mpiutil.gather_list(result_absorb, root = None)
        else:
            result_absorb = mpiutil.gather_list(result_absorb, root = None)
        if rank == 0:
            if self.critical_dis == True:
                with h5py.File(str(self.v)+'MHz_critical_dist.hdf5','w') as f:
                    f.create_dataset('critical_distance',data = result_absorb)
                    print ('critical_dis is saved')
            #if self.test == True:
            if False:
                with h5py.File('./' + str(self.emi_form)+str(self.v) + 'F2py_absorb.hdf5', 'w') as f:
                    f.create_dataset('F2py_absorb', data = result_absorb)
            else:
                result_absorb = np.array(result_absorb)
                absorb = result_absorb[:,1]
                I_E = self.I_E(self.v)
                result = []
                print ('in the beginning')
                for pix_number in range(unabsorb.size):
                    n = np.arange(6,0,-1)
                    for num in n:
                        if int(unabsorb.size - pix_number) == int(10**num):
                            print ('left number of pixel',unabsorb.size - pix_number) 
                    X = unabsorb[pix_number] - self.I_E(self.v)
                    l, b = hp.pix2ang(self.nside, pix_number, nest = False, lonlat = True)
                    tao = self.Fortran2Py_optical_deepth(l, b)
                    Y = absorb[pix_number] - I_E * np.exp(-tao[-1])
                    mean_exptao = Y / X
                    pix_value = diffuse_raw_[pix_number] * mean_exptao + I_E*np.exp(-tao[-1])
                    result.append([pix_number,pix_value])
                    #print 'pixel_number',pix_number
                with h5py.File('./' + str(self.emi_form)+str(self.v)+'MHz_sky_map_with_absorption.hdf5','w') as h:
                    h.create_dataset('data',data = np.array(result)[:,1])
                    if self.index_type == 'constant_index':
                        index = self.Beta_G_constant
                        index = np.ones(12*self.nside*self.nside) * index
                    if self.index_type == 'freq_dependent_index':
                        index = self.Beta_G_freq
                        index = np.ones(12*self.nside*self.nside) * index
                    if self.index_type == 'direction_dependent_index':
                        index = self.Beta_G
                    h.create_dataset('spectral_index',data = index)
                    h.create_dataset('smooth_absorb',data = absorb)
                    print ('end, good job!, you are the best')

            return result
