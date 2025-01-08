# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 12:35:17 2025

@author: ADMIN
"""

## Visualizing the power components
## Display and save the RGB maps 

# =============================================================================
# The min-max value of each RGB channel is set differently depending upon the
# distinct RCS produced by these targets. For example, the "dihedral-like" (Pdl) 
# power is passed through the red color gun. Since "dihedral-like" produces strong
# backscatter, their r_max value is kept high. Similarly, "surface-like" targets 
# generally produce low backscatter; hence, b_max is kept relatively low compared 
# to other max ranges for red and green. 
# Users may change values accordingly to visualize the dual-pol scattering power components better.
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import os

def RGB_func(Pd1, Pu1, Ps1, fdir, method_p, p_type1):
    
    if method_p == "dcmp":
        
        Pd1 = np.nan_to_num(Pd1, nan=0)
        Pu1 = np.nan_to_num(Pu1, nan=0)
        Ps1 = np.nan_to_num(Ps1, nan=0)
        
        r_min, r_max = 0, 0.25
        g_min, g_max = 0, 0.15
        b_min, b_max = 0, 0.10
    
        # Function to normalize each channel to the [0, 1] range based on its specific range
        def normalize_to_1(arr, min_val, max_val):
            return (arr - min_val) / (max_val - min_val)

        # Normalize each channel to the [0, 1] range based on its own range
        R_normalized = normalize_to_1(Pd1, r_min, r_max)
        G_normalized = normalize_to_1(Pu1, g_min, g_max)
        B_normalized = normalize_to_1(Ps1, b_min, b_max)

        # Stack the normalized channels to form the RGB image
        rgb_image = np.stack((R_normalized, G_normalized, B_normalized), axis=-1)
        rgb_image = np.clip(rgb_image, 0, 1)

        # Display the resulting image with proper scaling
        plt.imshow(rgb_image)
        plt.axis('off')  # Hide axes
        
        os.chdir(fdir)
        plt.savefig('DP_Powers_Dcmp_RGB_'+p_type1, dpi=300, bbox_inches='tight')
    

    else:
        
        Pd1 = np.nan_to_num(Pd1, nan=0)
        Pu1 = np.nan_to_num(Pu1, nan=0)
        Ps1 = np.nan_to_num(Ps1, nan=0)
        
        r_min, r_max = 0, 0.15
        g_min, g_max = 0, 0.10
        b_min, b_max = 0, 0.06
    
        # Function to normalize each channel to the [0, 1] range based on its specific range
        def normalize_to_1(arr, min_val, max_val):
            return (arr - min_val) / (max_val - min_val)

        # Normalize each channel to the [0, 1] range based on its own range
        R_normalized = normalize_to_1(Pd1, r_min, r_max)
        G_normalized = normalize_to_1(Pu1, g_min, g_max)
        B_normalized = normalize_to_1(Ps1, b_min, b_max)

        # Stack the normalized channels to form the RGB image
        rgb_image = np.stack((R_normalized, G_normalized, B_normalized), axis=-1)
        rgb_image = np.clip(rgb_image, 0, 1)

        # Display the resulting image with proper scaling
        plt.imshow(rgb_image)
        plt.axis('off')  # Hide axes
        
        os.chdir(fdir)
        plt.savefig('DP_Powers_fact_RGB_'+p_type1, dpi=300, bbox_inches='tight')
        
        
        