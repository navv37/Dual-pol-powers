# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 10:44:17 2024

Description: This codes computes the scattering power componennts for 
SLC (intensity and phase) dual-pol SAR data (decomposition based approach)

@author: Mr. Abhinav Verma and Prof. Avik Bhattacharya,
MRSLab, CSRE, IIT Bombay

Related publication: Abhinav Verma, Avik Bhattacharya, Subhadip Dey, Armando Marino, and Paolo Gamba 2024 
"Target Characterization and Scattering Power Components From Dual-Pol Sentinel-1 SAR Data".
IEEE Transactions on Geoscience and Remote Sensing, vol. 62, pp. 1-19, 2024, Art no. 5222619.
DOI: 10.1109/TGRS.2024.3460476
 
List of inputs:
1. A processed C2 matrix
2. A slope file (tool avilable in SNAP software)
Note 1 and 2 should be present inside a same folder
    
"""

import spectral
import numpy as np
import os
import warnings
from scipy.signal import convolve2d

sdir = os.getcwd()

n = input("Please enter (or paste) the complete file directory: ")
os.chdir(n)

method = "dcmp"
p_type = "SLC"

warnings.filterwarnings("ignore")

##### Reading the C2 matrix
C11_image = spectral.open_image('C11.bin.hdr')
C11 = C11_image.load()
C11 = C11.reshape(C11.shape[0], -1)

C12_image_r = spectral.open_image('C12_real.bin.hdr')
C12_r = C12_image_r.load()
C12_r = C12_r.reshape(C12_r.shape[0], -1)

C12_image_i = spectral.open_image('C12_imag.bin.hdr')
C12_i = C12_image_i.load()
C12_i = C12_i.reshape(C12_i.shape[0], -1)

C22_image = spectral.open_image('C22.bin.hdr')
C22 = C22_image.load()
C22 = C22.reshape(C22.shape[0], -1)

## Comment if slope mask is not required
Slope_image = spectral.open_image('slope.bin.hdr')
Slope = Slope_image.load()
Slope = Slope.reshape(Slope.shape[0], -1)

##### Averaging the values over a square kernel
ws = int(input("Please enter window size (recommended values are 5 or 7): "))
kernel = np.ones((ws,ws),np.float32)/(ws*ws);
pad = int((ws-1)/2);


C11_av = convolve2d(C11, kernel, mode = 'valid');
C11_av = np.pad(C11_av, pad_width=pad, mode='constant', constant_values=0)

C11_av_db = 10*np.log10(C11_av);

C12_r_av = convolve2d(C12_r, kernel, mode = 'valid');
C12_r_av = np.pad(C12_r_av, pad_width=pad, mode='constant', constant_values=0)

C12_i_av = convolve2d(C12_i, kernel, mode = 'valid');
C12_i_av = np.pad(C12_i_av, pad_width=pad, mode='constant', constant_values=0)

# C12_av = C12_r_av + 1j*C12_i_av;
# C21_av = np.conj(C12_av);

C22_av = convolve2d(C22, kernel, mode = 'valid');
C22_av = np.pad(C22_av, pad_width=pad, mode='constant', constant_values=0)

Slope_av = convolve2d(Slope, kernel, mode = 'valid');
Slope_av = np.pad(Slope_av, pad_width=pad, mode='constant', constant_values=0)

Slope_msk = np.where(Slope_av >= 12, 0, 1);

##### Calculating Stokes vector elements

s0 = C11_av + C22_av;
s1 = C11_av - C22_av;
s2 = 2*C12_r_av;
s3 = 2*C12_i_av;

##### Calculate Entropy
## Here eigen values are calculated using Stokes vector elements

tpp = np.sqrt(np.square(s1) + np.square(s2) + np.square(s3));

lmbd1 = (s0 + tpp)/2;
lmbd2 = (s0 - tpp)/2;

prob1 = lmbd1/(lmbd1 + lmbd2);
prob2 = lmbd2/(lmbd1 + lmbd2);

ent = -prob1*np.log2(prob1) - prob2*np.log2(prob2);

dop = (lmbd1 - lmbd2)/(lmbd1 + lmbd2)
beta = lmbd1/(lmbd1 + lmbd2)

##### Taking abs of Stokes vector elements
s0 = np.abs(s0);
s1 = np.abs(s1);
s2 = np.abs(s2);
s3 = np.abs(s3);

##### Slope mask the Stokes vector elements
s1_mask = np.multiply(s1,Slope_msk);
s2_mask = np.multiply(s2,Slope_msk);
s3_mask = np.multiply(s3,Slope_msk);

##### Normalizing Stokes vector elements

def S_norm(S_array):
    S_5 = np.percentile(S_array, 5)
    S_95 = np.percentile(S_array, 95)
    S_cln = np.where(S_array > S_95, S_95, S_array)
    S_cln = np.where(S_cln < S_5, S_5, S_cln)
    S_cln_max = np.max(S_cln)
    S_norm_array = np.divide(S_cln,S_cln_max) 
    
    return S_norm_array

s1_norm = S_norm(s1_mask)
s2_norm = S_norm(s2_mask)
s3_norm = S_norm(s3_mask)

s1_s_norm = S_norm(s1) #This is S1 normalzied for DpRSI, does not include slope mask

##### Power Calculation

dprbi = np.sqrt(np.square(s1_norm) + np.square(s2_norm) + np.square(s3_norm))/np.sqrt(3);

dprsi_con1 = (1 - ent)*np.sqrt(1 - np.square(s1_s_norm)); # For Valid pixels
dprsi_con2 = np.sqrt(1 - np.square(s1_s_norm)); # For Noise pixels 

NESZ = -16 ## For Sentinel-1
dprsi = np.where(C11_av_db > NESZ, dprsi_con1, dprsi_con2) 

shp = np.shape(dprbi)


alpha1 = np.arctan2(dprbi, 1 - dprbi)
alpha1 = np.degrees(alpha1)

alpha2 = np.arctan2(1-dprsi, dprsi)
alpha2 = np.degrees(alpha2)

alpha_dp = (alpha1 + alpha2)/2; #Dual-pol target characteristic parameter proposed in Verma et al. 2024

## Alpha as geomteric factor
alpha_dp_rad = np.radians(2*alpha_dp) 

cos_a = np.cos(alpha_dp_rad)

## Power components for valid pixels (VV > NESZ)
Pu_v = (1 - dop)*s0 
Pd_v = (1/2)*dop*s0*(1 - cos_a)
Ps_v = (1/2)*dop*s0*(1 + cos_a)

## Power components for noise pixels (VV < NESZ)
Pu_n = (1 - beta)*s0 
Pd_n = (1/2)*beta*s0*(1 - cos_a)
Ps_n = (1/2)*beta*s0*(1 + cos_a)

## Dual-pol scattering power
Pu = np.where(C11_av_db > NESZ, Pu_v, Pu_n) # Unpolized power
Pd = np.where(C11_av_db > NESZ, Pd_v, Pd_n) # "Dihedral-like" power
Ps = np.where(C11_av_db > NESZ, Ps_v, Ps_n) # "Surface-like" power

## Saving the output

dprbi = dprbi.astype(np.float32)
dprsi = dprsi.astype(np.float32)
Pd = Pd.astype(np.float32)
Ps = Ps.astype(np.float32)
Pu = Pu.astype(np.float32)

metadata = C11_image.metadata

metadata.update({"band names":"DpRBI_slc"})
spectral.envi.save_image('DpRBI_slc.hdr', dprbi, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"DpRSI_slc"})
spectral.envi.save_image('DpRSI_slc.hdr', dprsi, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"Alpha_dp_slc"})
spectral.envi.save_image('Alpha_dp_slc.hdr', alpha_dp, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"Pdl_dcmp_slc"})
spectral.envi.save_image('Pdl_dcmp_slc.hdr', Pd, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"Psl_dcmp_slc"})
spectral.envi.save_image('Psl_dcmp_slc.hdr', Ps, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"Pu_dcmp_slc"})
spectral.envi.save_image('Pu_dcmp_slc.hdr', Pu, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)


## RGB display of the dual-pol powers
os.chdir(sdir)

import RGB_display_dual_pol_powers
RGB_display_dual_pol_powers.RGB_func(Pd, Pu, Ps, n, method, p_type)





