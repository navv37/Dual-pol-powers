# Dual-pol scattering power components 
Author - Mr. Abhinav Verma and Prof. Avik Bhattacharya

MRSLab, Indian Institute of Technology Bombay, India

## Method-I: Factorization based approach

**Publication:** Abhinav Verma, Avik Bhattacharya, Subhadip Dey, Carlos López-Martínez, and Paolo Gamba, 
"Scattering power components from dual-pol Sentinel-1 SLC and GRD SAR data".
*ISPRS Journal of Photogrammetry and Remote Sensing*, vol.212, pp. 289-305, 2023.
DOI: https://doi.org/10.1016/j.isprsjprs.2024.05.010

**To get the dual-pol scattering power components using Method-I:**
1. Use Python code "Dual_pol_powers_factorization_slc.py" for SLC (intensity and phase information) SAR data
2. Use Python code "Dual_pol_powers_factorization_grd.py" for GRD (only intensity, no phase information) SAR data

## Method-II: Decomposition based approach

**Publication:** Abhinav Verma, Avik Bhattacharya, Subhadip Dey, Armando Marino, and Paolo Gamba, 
"Target Characterization and Scattering Power Components From Dual-Pol Sentinel-1 SAR Data".
*IEEE Transactions on Geoscience and Remote Sensing*, vol. 62, pp. 1-19, 2024, Art no. 5222619.
DOI: https://doi.org/10.1109/TGRS.2024.3460476

**To get the dual-pol scattering power components using Method-II:**
1. Use Python code "Dual_pol_powers_decompostion_slc.py" for SLC (intensity and phase information) SAR data
2. Use Python code "Dual_pol_powers_decompostion_grd.py" for GRD (only intensity, no phase information) SAR data

### Note: Please keep "RGB_display_dual_pol_powers.py" with Method-I and Method-II Python files to generate RGB images of the powers.

## Requirements

**Python libraries required:**

1. spectral
2. numpy
3. scipy
4. matplotlib

**List of inputs to run the code:**

 A. For SLC product: Elements of the processed C2 matrix alongside a slope file (tool available in SNAP software)
  1. C11.bin.hdr
  2. C12_real.bin.hdr
  3. C12_imag.bin.hrd
  4. C22.bin.hdr
  5. slope.bin.hdr

 B. For GRD product: Processed sigma0_VV and sigma0_VH file
  1. Sigma0_VV.bin.hdr
  2. Sigma0_VH.bin.hdr
  3. slope.bin.hdr

**Standard processing steps to process dual-pol SAR data in SNAP:**

 A. For SLC product
  1. TopSAR split 
  2. Apply "orbit file"
  3. Radiometric calibration  
  4. TopSAR Deburst
  5. Polarimteric matrix [C2]
  6. Multi-looking
  7. Polarimteric speckle filter (optional)
  8. Compute slope (Raster -> DEM Tools -> Compute Slope and Aspect)
  9. Export as "PolSARpro" format (File -> Exports -> SAR Formats -> PolSARPro)

 B. For GRD product
  1. Apply "orbit file"
  2. Radiometric calibration
  3. Speckle filter
  4. Compute slope (Raster -> DEM Tools -> Compute Slope and Aspect)
  5. Export as "PolSARpro" format (File -> Exports -> SAR Formats -> PolSARPro)
