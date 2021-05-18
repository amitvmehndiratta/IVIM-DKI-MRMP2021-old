# IVIM-DKI analysis with total variation penalty funtion in prostate cancer
---------------------------------------------------------------------------------

|<https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021>|
--> A package providing MATLAB programming tools for IVIM-DKI analysis with total
variation (TV) penalty function.

---------------------------------------------------------------------------------
**REFERENCES:** [1] Kayal, E. B. et al. (2017). Quantitative analysis of intravoxel 
incoherent motion (IVIM) diffusion MRI using total variation and Huber penalty function. 
Medical physics, 44(11), 5849-5858.
[2] Malagi, A. V. et al. (2021). IVIM-DKI differentiation between prostate cancer 
and benign prostatic hyperplasia: comparison of 1.5T vs. 3T MRI. 
Magnetic Resonance Materials in Physics, Biology and Medicine.
---------------------------------------------------------------------------------
**AUTHOR:** Dr. Amit Mehndiratta, Indian Institute of Technology Delhi, India 
E-mail: <amehndiratta@cbme.iitd.ac.in>, <amit.mehndiratta@keble.oxon.org>.
---------------------------------------------------------------------------------
**Disclaimer: This project can be used only for research purposes. Authors are 
not liable for any clinical use of it, authors could not be held responsible.**
---------------------------------------------------------------------------------
- Version 1.0: May 2021
---------------------------------------------------------------------------------

- If you have any queries or suggestions about this package, 
    please do not hesitate to contact us.


### This package consists of two folders:

[1] **'IVIM_DKI_data'**: The folder consists of IVIM-DKI 4D data with tumor, BPH, 
and healthy PZ masks. All the images and ROIs are in compressed nii format and to 
read the images either 'niftiread' function (Matlab R2017b and later) or 
'load_untouch_niigz' function from
https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
'main_program.m': Main program which executes all funtions and codes together.

2. **'IVIM_DKI_functions'**: MATLAB codes to compute IVIM-DKI parameter maps obtained from 
traditional IVIM-DKI and novel IVIM-DKI model with TV. 
'AIC.m': Calculates Akaike information criterion (AIC) and AIC corrected.
'aic_aicccal.m': Calculates AIC/AICc for tumor, BPH and healthy PZ masks.
'allivimdki.m': Function which contains IVIM-DKI model equation.
'hybrid_TVmodel.m': Executes IVIM-DKI model with TV code.
'hybridmodel.m': Executes Monoexponential and IVIM-DKI model code. 
'monoexplog.m': Function which contains monoexponential model equation.
'tv3d.m':Function which contains 3D TV.
Please see ref. [1] for more details on the implementation of TV function.

**Note: IVIM_DKI_functions folder must be setpath before execution of 'main_program.m'.**