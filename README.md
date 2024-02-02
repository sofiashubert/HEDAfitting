# HEDAfitting

MATLAB code to fit the extinction spectra of plasmonic colloidal semiconductor nanocrystals(NCs).

GitHub repository, https://github.com/sofiashubert/HEDAfitting

This code is open source and available for use by anyone. When publishing please
use the following citations:
Gibbs, S.L. et al. JPCC (2020) 
Shubert-Zuleta, S.A. et al. (2024)

This code was originally developed by Dr. Stephen Gibbs and Dr. Corey Staller to fit the plasmon peak of Sn:In2O3 NCs. It was edited for clarity and efficiency by Dr. Zachary Sherman and later adapted for nonparabolic band structure materials by Sofia Shubert-Zuleta. Questions regarding the fitting can be directed to sofiashubert@gmail.com

The MATLAB code is meant to fit multiple spectra at once, in this example it fits spectra all within one reducing titration (same NCs in each spectra but reduced to different extents) but the code could also be used to fit the spectra of many unrelated spectra. The input file (test.txt) includes the spectral data with the first column being wavenumbers and the subsequent columns being the extinction data. The number of columns with spectral data must match the length of the input arrays “CoCp2_vol”, “a_mean”, and “a_dev”. The code outputs the fit results into an excel file.

The mean radius (a_mean), standard deviation in the radius (a_dev), and dispersion volume fractions (eta) must be input from experimental data in the arrays noted in parentheses. The following variables are constants that must be changed depending on the NC material, solvent properties, and experimental parameters. 

m_eff_init = effective mass of the as-synthesized NCs
m_eff_CBmin = effective mass at the conduction band minimum (only used for nonparabolic materials)
C_nonparabolic = nonparabolicity factor for the material. 
A_surf = the material dependent proportionality factor to describe surface damping (typically 0-2)
Path_length = path length of the cuvette used for the measurement. Units are in meters.
eps_inf = high frequency dielectric constant of the NC material
eps_f = solvent dielectric constant

The above describes fitting multiple spectra for a nonparabolic conduction band material where the effective mass is changing between spectra but the initial effective mass of the first spectrum is a known value. If the NC material has a parabolic conduction band (like Sn:In2O3 NCs). Line 129 can be commented out and replaced with the contents in line 128, which populates the effect mass matrix with the known effective mass of the NCs. 
