
# coding: utf-8
#!/usr/bin/env python
import matplotlib
#get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import numpy as np;import healpy as hp;import h5py;
import scipy.constants as Cons
K = Cons.k
C = Cons.c
from scipy.special import erf
import sys
import os

class direction_dependent_index(object):
    def __init__(self,v_file_dir=None):
        _path = os.path.dirname(os.path.abspath(__file__))
        _path = os.path.split(_path)[0]
        _path = os.path.split(_path)[0]
        self.file_dir1 = _path +'/ULSA/spectral_index_fitting'
        self.v_file_dir = v_file_dir
        self.file_dir = _path +'/obs_sky_data'
    def unit(self,v):
        return np.square(v*1e6) * 2*K/C**2
    
    def change_coord(self, m, coord):
        """ Change coordinates of a HEALPIX map

        Parameters
        ----------
        m : map or array of maps
          map(s) to be rotated
        coord : sequence of two character
          First character is the coordinate system of m, second character
          is the coordinate system of the output map. As in HEALPIX, allowed
          coordinate systems are 'G' (galactic), 'E' (ecliptic) or 'C' (equatorial)

        Example
        -------
        The following rotate m from galactic to equatorial coordinates.
        Notice that m can contain both temperature and polarization.
        >>>> change_coord(m, ['G', 'C'])
        """
        # Basic HEALPix parameters
        npix = m.shape[-1]
        nside = hp.npix2nside(npix)
        ang = hp.pix2ang(nside, np.arange(npix))
    
        # Select the coordinate transformation
        rot = hp.Rotator(coord=reversed(coord))

        # Convert the coordinates
        new_ang = rot(*ang)
        new_pix = hp.ang2pix(nside, *new_ang)

        return m[..., new_pix]

    def minus_free_free(self,m,v,downgrade_to):
        #in galactic coordinate
        with h5py.File(self.file_dir + '/Free_free_emission/T_ff_nu_'+str(v)+'MHz.hdf5','r') as f:
            free_free_G = f['free_free'][:]
        #free_free_G = hp.ud_grade(free_free_G,downgrade_to)
        
        #smothing the data in order to have same resolution.
        free_free_G = self.smooth(free_free_G)
        m = self.change_coord(m,["C","G"])
        m = self.masked_smoothing(m)
        m = self.change_coord(m,["G","C"])

        free_free_G = self.change_coord(free_free_G,['G','C']) 
        result = m - free_free_G
        #setting the value less than 0.1k to nan value.
        result[np.where(result<0.1)] = np.nan
        #print ('np.where(result<0.1',np.where(result<0.1)[0].shape)
        return result

    
    def IndexToDeclRa(self,index,downgrade_to):
        theta,phi=hp.pix2ang(downgrade_to,index)
        return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)

    def read_file(self,downgrade_to= 256, resolution = 5):

        with h5py.File(self.file_dir + '/Guzman/wlb_45MHz.hdf5','r') as h:
            hpmap_45_old = h['hpmap'][:]
        hpmap_45_old = self.minus_free_free(hpmap_45_old,45,downgrade_to)
        hpmap_408 = hp.read_map(self.file_dir + '/Haslam/haslam408_dsds_Remazeilles2014.fits')
        #downgrade to 256 may be, and change the coordinate from galaxy to equatorial
        hpmap_408 = hp.ud_grade(hpmap_408,downgrade_to)
        hpmap_408 = self.change_coord(hpmap_408,['G', 'C'])
        if True:
            hpmap_408 = self.minus_free_free(hpmap_408,408,downgrade_to)
        
        #####
        #the data from LWA    
        hpmap_35 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-35.fits')
        hpmap_35 = self.minus_free_free(hpmap_35,35,downgrade_to)
        hpmap_38 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-38.fits')
        hpmap_38 = self.minus_free_free(hpmap_38,38,downgrade_to)
        hpmap_40 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-40.fits')
        hpmap_40 = self.minus_free_free(hpmap_40,40,downgrade_to)
        hpmap_45 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-45.fits')
        hpmap_45 = self.minus_free_free(hpmap_45,45,downgrade_to)
        hpmap_50 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-50.fits')
        hpmap_50 = self.minus_free_free(hpmap_50,50,downgrade_to)
        hpmap_60 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-60.fits')
        hpmap_60 = self.minus_free_free(hpmap_60,60,downgrade_to)
        hpmap_70 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-70.fits')
        hpmap_70 = self.minus_free_free(hpmap_70,70,downgrade_to)
        hpmap_74 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-74.fits')
        hpmap_74 = self.minus_free_free(hpmap_74,74,downgrade_to)
        hpmap_80 = hp.read_map(self.file_dir + '/LWA/healpix-all-sky-rav-rsclean-map-80.fits')
        hpmap_80 = self.minus_free_free(hpmap_80,80,downgrade_to)
        #print ('np.isnan(hpmap_35)',np.isnan(hpmap_35)) 
        return hpmap_408,hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80
 

       

    
    def masked_smoothing(self,U, rad=5.0):     
        V=U.copy()
        V[U!=U]=0
        VV=hp.smoothing(V, fwhm=np.radians(rad))    
        W=0*U.copy()+1
        W[U!=U]=0
        WW=hp.smoothing(W, fwhm=np.radians(rad))    
        return VV/WW
    

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
    
    def smooth(self,map_in, fwhm = 5.0):
        fwhm = np.radians(fwhm)
        return hp.smoothing(map_in,fwhm)


    
    def smoothing_data(self,downgrade_to):
        hpmap_408, hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        
        
        
        Dict = {'hpmap_408':hpmap_408,'hpmap_45_old':hpmap_45_old,'hpmap_35':hpmap_35,'hpmap_38':hpmap_38,'hpmap_40':hpmap_40,'hpmap_45':hpmap_45,'hpmap_50':hpmap_50,'hpmap_60':hpmap_60,'hpmap_70':hpmap_70,'hpmap_74':hpmap_74,'hpmap_80':hpmap_80}
        for key,X in Dict.items():
            X = self.change_coord(X,['C','G'])
            nans, x= self.nan_helper(X)
            X[nans]= np.interp(x(nans), x(~nans), X[~nans])
            #smooth step have been done in minus_free_free
            #X = self.smooth(X)

            X = hp.ud_grade(X,downgrade_to)
            if key == 'hpmap_35':
                Mask_missing_region_lwa = nans.copy()
                Mask_missing_region_lwa = hp.ud_grade(Mask_missing_region_lwa,downgrade_to)
                
            Dict[key] = X
            
        #the output map with coordinate of galaxy
        return Dict,Mask_missing_region_lwa
        
    def I_E(self,v):
        #result = 24.4 *(v*1e-3/0.31)**-2.58
        result = 1.2*(v*1e-3/1.0)**-2.58
        #print (result)
        return result
    def func1(self,beta, x1, x2):
        return (x2-self.I_E(408.)) * (x1/408.)**beta

    def error1(self,beta, x1, x2, y):
        #print 'error1','x1,x2,y.shape',x1.shape,x2.shape,y.shape
        return (self.func1(beta,x1,x2) - (y-self.I_E(x1)))/(y) 
        #return (func1(beta,x1,x2) - (y-I_E(x1)))


    def calculate_index(self,downgrade_to):
        
        Dict,Mask_missing_region_lwa = self.smoothing_data(downgrade_to)
        hpmap_408, hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = Dict['hpmap_408'],Dict['hpmap_45_old'],Dict['hpmap_35'],Dict['hpmap_38'],Dict['hpmap_40'],Dict['hpmap_45'],Dict['hpmap_50'],Dict['hpmap_60'],Dict['hpmap_70'],Dict['hpmap_74'],Dict['hpmap_80']
        
        #hpmap_408, hpmap_45_old,hpmap_35,hpmap_38,hpmap_40,hpmap_45,hpmap_50,hpmap_60,hpmap_70,hpmap_74,hpmap_80 = self.read_file()
        hpmap_40 = self.change_coord(hpmap_40,['G','C'])
        hpmap_45 = self.change_coord(hpmap_45,['G','C'])
        hpmap_50 = self.change_coord(hpmap_50,['G','C'])
        hpmap_60 = self.change_coord(hpmap_60,['G','C'])
        hpmap_70 = self.change_coord(hpmap_70,['G','C'])
        pix_number_45 = []

        for dec in np.arange(-40,0,0.1):
            for ra in np.arange(-100,-30,0.1):
                pix_num = self.DeclRaToIndex(dec,ra,downgrade_to)
                pix_number_45.append(pix_num)
        
        # at this position the pixel value - I_E(v) nearly zero or less zero, so we dont need the data when calcule index in that LOS, by setting the value to zero, so it can be tell apart by the mask_condition in line 226
        hpmap_40[pix_number_45] = 0. 
        hpmap_45[pix_number_45] = 0. 
        hpmap_50[pix_number_45] = 0. 
        hpmap_60[pix_number_45] = 0.
        hpmap_70[pix_number_45] = 0. 
        hpmap_40 = self.change_coord(hpmap_40,['C','G'])
        hpmap_45 = self.change_coord(hpmap_45,['C','G'])
        hpmap_50 = self.change_coord(hpmap_50,['C','G'])
        hpmap_60 = self.change_coord(hpmap_60,['C','G'])
        hpmap_70 = self.change_coord(hpmap_70,['C','G'])
        
        #freq = np.array([22,35,38,40,45,50,60,70,74,80])
        freq = np.array([35,38,40,45,50,60,70,74,80])
        #extra data part
        if self.v_file_dir != None:
            result = {}

            for v,dir_ in self.v_file_dir.items():
                f = h5py.File(self.file_dir + dir_,'r')
                data = f['data'][:]
                X = data.copy()
                X = self.change_coord(X,['C','G'])
                nans, x= self.nan_helper(X)
                X[nans]= np.interp(x(nans), x(~nans), X[~nans])
                X = self.smooth(X)
                #so the nan value region will be tell apart by >0 rule
                X[nans] = 0.
                X = hp.ud_grade(X,downgrade_to)
                f.close()
                result[v] = X

            Freq = []
            Value = []
            for key,value in result.items():
                #print (eval(key))
                Freq.append(key)
                Value.append(value)
            freq = list(freq) + Freq
            freq = np.array(freq)

        spectral_index_lwa_and_408 = []
        for i in range(12*downgrade_to**2):
            mask_condition = np.array([hpmap_35[i]-2.725-self.I_E(35), hpmap_38[i]-2.725-self.I_E(38), hpmap_40[i]-2.725-self.I_E(40), hpmap_45[i]-2.725-self.I_E(45), hpmap_50[i]-2.725-self.I_E(50), hpmap_60[i]-2.725-self.I_E(60), hpmap_70[i]-2.725-self.I_E(70), hpmap_74[i]-2.725-self.I_E(74), hpmap_80[i]-2.725-self.I_E(80)])
            if self.v_file_dir != None:
                result_ = []
                for j in range(len(Value)):
                    result_.append(Value[j][i]-2.725-self.I_E(Freq[j]))
                mask_condition = list(mask_condition) + result_
                mask_condition = np.array(mask_condition)

            mask_ = np.where(mask_condition>0)[0]
            value_freq = np.array([hpmap_35[i]-2.725, hpmap_38[i]-2.725, hpmap_40[i]-2.725, hpmap_45[i]-2.725, hpmap_50[i]-2.725, hpmap_60[i]-2.725, hpmap_70[i]-2.725, hpmap_74[i]-2.725, hpmap_80[i]-2.725])
            if self.v_file_dir != None:
                append_value_freq = []
                for j in range(len(Value)):
                    append_value_freq.append(Value[j][i]-2.725)

                value_freq = list(value_freq) + append_value_freq
                value_freq = np.array(value_freq)

            value_freq = value_freq[mask_]
            
            
            value_408 = np.ones_like(value_freq) * hpmap_408[i]
            x1 = freq[mask_].copy()
            x2 = value_408.copy()
            y = value_freq.copy()
            
            beta = [-2.6]
            if int(value_freq.shape[0])==0:
                Para_constant = np.nan
                #print ('the ith pixel no available data',i)
            else:  
                Para_constant=leastsq(self.error1,beta,args=(x1,x2,y))[0][0]
            if i % 100 == 0:
                #print ('Para_constant',Para_constant,'left_number',12*256**2 - i)
            spectral_index_lwa_and_408.append(Para_constant)
        spectral_index_lwa_and_408 = np.array(spectral_index_lwa_and_408).reshape(-1)
        index_45_old_and_408,Mask = self.index_between_45_and_408(hpmap_45_old,hpmap_408)
        return spectral_index_lwa_and_408,index_45_old_and_408,Mask,Mask_missing_region_lwa
    
    def IndexToDeclRa(self,index,downgrade_to):
        theta,phi=hp.pix2ang(downgrade_to,index)
        return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)
    
    def DeclRaToIndex(self,decl,RA,downgrade_to):
        return hp.pixelfunc.ang2pix(downgrade_to,np.radians(-decl+90.),np.radians(360.-RA))
    
    def combined_index(self,downgrade_to):
        spectral_index_lwa_and_408,index_45_old_and_408,Mask,Mask_missing_region_lwa = self.calculate_index(downgrade_to)
        
        
        hpmap_45_old = self.change_coord(index_45_old_and_408,['G','C']).copy()
        hpmap_45 = self.change_coord(spectral_index_lwa_and_408,['G','C']).copy()
        
        new_map = []
        #LWA = -40; Guzman = 67
        Dec_0 = -15;LWA_bottom_limit = -40
        A =2. / ((Dec_0 - LWA_bottom_limit)/Dec_0)
        #for pix_number in range(12*downgrade_to**2):
        for pix_number in range(12*256**2):
            Dec,Ra = self.IndexToDeclRa(pix_number,downgrade_to) 
            pix_value = 0.5*(1 + erf(A*(Dec-Dec_0)/Dec_0))*hpmap_45[pix_number] + 0.5*(1 - erf(A*(Dec-Dec_0)/Dec_0))*hpmap_45_old[pix_number]
            new_map.append(pix_value)
        new_map = np.array(new_map)
        new_map = self.change_coord(new_map,['C','G'])
        nans, x = self.nan_helper(new_map)
        new_map[nans]=np.interp(x(nans), x(~nans), new_map[~nans])

        with h5py.File(self.file_dir1 + '/spectral_index_map.hdf5','w') as f:
            f['spectral_index'] = new_map
            f.close()

        #f1 = h5py.File(self.file_dir + '/spectral_index_map.hdf5', 'r+') # open the file
        #data = f1['spectral_index']       # load the data
        #data[...] = new_map               # assign new values to data
        #f1.close()       
        return new_map    
        
    def index_between_45_and_408(self,map_45_old,map_408):
        #interpolate the pixel less than I_E(45)
        mask = np.where(map_45_old - 2.725-self.I_E(45)<=0)[0]
        map_45_old[mask] = np.nan
        nans, x = self.nan_helper(map_45_old)
        map_45_old[nans]=np.interp(x(nans), x(~nans), map_45_old[~nans])
        
        mask_45_less_I_E = mask.copy()
        
        index = np.log10((map_45_old-2.725-self.I_E(45))/(map_408-self.I_E(408)))/np.log10(45./408.)
        index = index.reshape(-1)
        return index,mask_45_less_I_E
        

