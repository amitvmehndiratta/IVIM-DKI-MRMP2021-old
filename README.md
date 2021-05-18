# IVIM-DKI analysis with total variation penalty funtion in prostate cancer
---------------------------------------------------------------------------------

|<https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021>|
--> A package providing MATLAB programming tools for IVIM-DKI analysis with total
variation (TV) penalty function.

---------------------------------------------------------------------------------
**REFERENCES:** 
1. Kayal, E. B. et al. (2017). Quantitative analysis of intravoxel 
incoherent motion (IVIM) diffusion MRI using total variation and Huber penalty function. 
Medical physics, 44(11), 5849-5858.

2. Malagi, A. V. et al. (2021). IVIM-DKI differentiation between prostate cancer 
and benign prostatic hyperplasia: comparison of 1.5T vs. 3T MRI. 
Magnetic Resonance Materials in Physics, Biology and Medicine.
---------------------------------------------------------------------------------
**AUTHOR:** Dr. Amit Mehndiratta, Indian Institute of Technology Delhi, India 
E-mail: <amehndiratta@cbme.iitd.ac.in>, <amit.mehndiratta@keble.oxon.org>.
- If you have any queries or suggestions about this package, 
    please do not hesitate to contact us.
---------------------------------------------------------------------------------
Disclaimer: This project can be used only for research purposes. Authors are not liable for any clinical use of it, authors could not be held responsible.
---------------------------------------------------------------------------------
- Version 1.0: May 2021
---------------------------------------------------------------------------------

### This package consists of two folders:

1. **[IVIM_DKI_data](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/tree/main/IVIM_DKI_data)**: The folders consist of IVIM-DKI 4D data acquried at [1.5T](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/tree/main/IVIM_DKI_data/IVIM_DKI_1_5T) and [3T](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/tree/main/IVIM_DKI_data/IVIM_DKI_3T) MRI with tumor, BPH, and healthy PZ ROI masks. All the images and ROIs are in compressed nii format and to read these images either use 'niftiread' function (Matlab R2017b and later) or [load_untouch_niigz](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).


2. **[IVIM_DKI_functions](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/tree/main/IVIM_DKI_functions)**: MATLAB codes to compute IVIM-DKI parameter maps obtained from 
traditional IVIM-DKI and novel IVIM-DKI model with TV. 

-'AIC.m': Calculates Akaike information criterion (AIC) and AIC corrected (AICc).
-'aic_aicccal.m': Calculates AIC/AICc for tumor, BPH and healthy PZ ROIs from ['AIC.m'](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/blob/main/IVIM_DKI_functions/AIC.m).
-'allivimdki.m': Function which contains IVIM-DKI model equation.
-'hybrid_TVmodel.m': Executes IVIM-DKI model with TV code using [non-linear least square fitting with iterative total variation (TV) penalty function to perform spatial homogeneity on IVIM-DKI parameter reconstruction](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1002/mp.12520).
-'hybridmodel.m': Executes Monoexponential and IVIM-DKI model code using ['monoexplog.m'](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/blob/main/IVIM_DKI_functions/monoexplog.m) and ['allivimdki.m'](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/blob/main/IVIM_DKI_functions/allivimdki.m). 
-'im2Y.m':Transforms functional image data (4D or 3D array) into data matrix VxM where V is the number of voxels and M is number of b-values.
-'monoexplog.m': Function which contains monoexponential model equation.
-'tv3d.m':Function which calculates 3D TV penalty.

Please see ref. [1] for more details on the implementation of TV function.

Acknowledgement: https://github.com/oliverchampion/IVIM_tools

**Note: Please execute 'main_program.m' in IVIM_DKI_functions folder.**
