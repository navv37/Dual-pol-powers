# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 18:06:19 2024

Description: This codes computes the scattering power componennts for 
GRD (only intensity no phase) dual-pol SAR data (factorization based approach)

@author: Mr. Abhinav Verma and Prof. Avik Bhattacharya,
MRSLab, CSRE, IIT Bombay

Related publication: Abhinav Verma, Avik Bhattacharya, Subhadip Dey, Carlos López-Martínez, and Paolo Gamba 2023 
"Scattering power components from dual-pol Sentinel-1 SLC and GRD SAR data".
 ISPRS Journal of Photogrammetry and Remote Sensing, vol.212, pp. 289-305.
 DOI: https://doi.org/10.1016/j.isprsjprs.2024.05.010
 
List of inputs:
1. Sigma0_VV and Sigma0_VH files
2. A slope file (tool avilable in SNAP software)
Note 1 and 2 should be present inside a same folder
    
"""

import spectral
import numpy as np
import os
import warnings
from scipy.signal import convolve2d

sdir = os.getcwd()

n = input("Please enter (or paste) the file directory: ")
os.chdir(n)

method = "fact"
p_type = "GRD"

warnings.filterwarnings("ignore")

##### Reading the C2 matrix
C11_image = spectral.open_image('Sigma0_VV.bin.hdr')
C11 = C11_image.load()
C11 = C11.reshape(C11.shape[0], -1)

C22_image = spectral.open_image('Sigma0_VH.bin.hdr')
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


C22_av = convolve2d(C22, kernel, mode = 'valid');
C22_av = np.pad(C22_av, pad_width=pad, mode='constant', constant_values=0)

Slope_av = convolve2d(Slope, kernel, mode = 'valid');
Slope_av = np.pad(Slope_av, pad_width=pad, mode='constant', constant_values=0)

Slope_msk = np.where(Slope_av >= 12, 0, 1);

##### Calculating Stokes vector elements

s0 = C11_av + C22_av;
s1 = C11_av - C22_av;

##### Calculate Entropy
## Here eigen values are calculated using Stokes vector elements

prob1 = C11_av/(C11_av + C22_av);
prob2 = C22_av/(C11_av + C22_av);

ent = -prob1*np.log2(prob1) - prob2*np.log2(prob2);

##### Taking abs of Stokes vector elements
s0 = np.abs(s0); # Already postive
s1 = np.abs(s1);

##### Slope mask the Stokes vector elements
C11_mask = np.multiply(C11_av,Slope_msk);
C22_mask = np.multiply(C22_av,Slope_msk);
s1_mask = np.multiply(s1,Slope_msk);

##### Normalizing Stokes vector elements

def S_norm(S_array):
    S_5 = np.nanpercentile(S_array, 5)
    S_95 = np.nanpercentile(S_array, 95)
    S_cln = np.where(S_array > S_95, S_95, S_array)
    S_cln = np.where(S_cln < S_5, S_5, S_cln)
    S_cln_max = np.nanmax(S_cln)
    S_norm_array = np.divide(S_cln,S_cln_max) 
    
    return S_norm_array

C11_norm = S_norm(C11_mask)
C22_norm = S_norm(C22_mask)
s1_norm = S_norm(s1_mask)

s1_s_norm = S_norm(s1) #This is S1 normalzied for DpRSI, does not include slope mask

##### Power Calculation

dprbi = np.sqrt(np.square(C11_norm) + np.square(C22_norm))/np.sqrt(2);
dprbi = dprbi*s1_norm;

dprsi_con1 = (1 - ent)*np.sqrt(1 - np.square(s1_s_norm)); # For Valid pixels
dprsi_con2 = np.sqrt(1 - np.square(s1_s_norm)); # For Noise pixels 

NESZ = -16 ## For Sentinel-1
dprsi = np.where(C11_av_db > NESZ, dprsi_con1, dprsi_con2) 

shp = np.shape(dprbi)

dprbi_flt = dprbi.flatten()
dprsi_flt = dprsi.flatten()
shp_flt = np.shape(dprbi_flt)


indices_vec = np.array([dprsi_flt, dprbi_flt]).transpose()
indices_vec_sort = np.array([[max(row), min(row)] for row in indices_vec])

y1 = indices_vec_sort[:,0] #First dominant
y2 = (1 - indices_vec_sort[:,0])*indices_vec_sort[:,1] #Second dominant

residue = 1 - (y1 + y2)

dprsi_dom = np.where(dprsi_flt > dprbi_flt)[0] #Keeps the tuple where dprsi is dominant
dprbi_dom = np.where(dprsi_flt < dprbi_flt)[0] #Keeps the tuple where dprbi is dominant

## Surface-like power component
Ps = np.zeros(shp_flt)
##dprsi_dom and dprbi_dom are not dprsi and dprbi values, they just indicate tuples (pixel) for which they are greater
Ps[dprsi_dom] = y1[dprsi_dom]; #In these tuples dprsi was dominant, hence taking y1 (first dominant)
Ps[dprbi_dom] = y2[dprbi_dom]; #In these tuples dprbi was dominant, hence taking y2 (Second domiant)

Ps = Ps.reshape(shp[0],shp[1])
Ps = np.multiply(s0,Ps)

## Dihedral-like power component
Pd = np.zeros(shp_flt)
Pd[dprbi_dom] = y1[dprbi_dom];
Pd[dprsi_dom] = y2[dprsi_dom];

Pd = Pd.reshape(shp[0],shp[1])
Pd = np.multiply(s0,Pd)

## Residue (diffused) power component

Pr = residue.reshape(shp[0],shp[1])
Pr = np.multiply(s0,Pr)

## Saving the output

dprbi = dprbi.astype(np.float32)
dprsi = dprsi.astype(np.float32)
Pd = Pd.astype(np.float32)
Ps = Ps.astype(np.float32)
Pr = Pr.astype(np.float32)

metadata = C11_image.metadata

metadata.update({"band names":"DpRBI_grd"})
spectral.envi.save_image('DpRBI_grd.hdr', dprbi, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"DpRSI_grd"})
spectral.envi.save_image('DpRSI_grd.hdr', dprsi, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"Pdl_fact_grd"})
spectral.envi.save_image('Pdl_fact_grd.hdr', Pd, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"Psl_fact_grd"})
spectral.envi.save_image('Psl_fact_grd.hdr', Ps, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

metadata.update({"band names":"Pr_fact_grd"})
spectral.envi.save_image('Pr_fact_grd.hdr', Pr, metadata=metadata, interleave = 'bsq', ext='.bin', force=True)

## RGB display of the dual-pol powers
os.chdir(sdir)

import RGB_display_dual_pol_powers
RGB_display_dual_pol_powers.RGB_func(Pd, Pr, Ps, n, method, p_type)



