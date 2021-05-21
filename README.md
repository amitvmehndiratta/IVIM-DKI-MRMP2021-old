# IVIM-DKI analysis with total variation penalty funtion in prostate cancer
---------------------------------------------------------------------------------

|<https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021>|
--> A package providing MATLAB programming tools for IVIM-DKI analysis with total
variation (TV) penalty function that provides an adaptive spatial homogeneity by 
removing spurious values that may arise during parametric reconstruction of the image.

---------------------------------------------------------------------------------
**REFERENCES:** 

[1](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1002/mp.12520). Kayal, E. B. et al. (2017). Quantitative analysis of intravoxel 
incoherent motion (IVIM) diffusion MRI using total variation and Huber penalty function. 
Medical physics, 44(11), 5849-5858.

[2]. Malagi, A. V. et al. (2021). IVIM-DKI differentiation between prostate cancer 
and benign prostatic hyperplasia: comparison of 1.5T vs. 3T MRI. 
Magnetic Resonance Materials in Physics, Biology and Medicine.

---------------------------------------------------------------------------------
If you have any queries or suggestions about this package, kindly contact 
Dr. Amit Mehndiratta, Indian Institute of Technology Delhi, India. 

E-mail: <amehndiratta@cbme.iitd.ac.in>, <amit.mehndiratta@keble.oxon.org>.

---------------------------------------------------------------------------------
Disclaimer: This project can be used only for research purposes. Authors are not liable for any clinical use of it, authors could not be held responsible.

---------------------------------------------------------------------------------
- Version 1.0: May 2021
---------------------------------------------------------------------------------

### This package consists of two folders:

1. **[IVIM_DKI_data](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/tree/main/IVIM_DKI_data)**: The folders consist of IVIM-DKI 4D data acquried at [1.5T](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/tree/main/IVIM_DKI_data/IVIM_DKI_1_5T) and [3T](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/tree/main/IVIM_DKI_data/IVIM_DKI_3T) MRI with tumor, BPH, and healthy PZ ROI masks. All the images and ROIs are in compressed nii format and to read these images either use 'niftiread' function (Matlab R2017b and later) or [load_untouch_niigz](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).


2. **[IVIM_DKI_functions](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/tree/main/IVIM_DKI_functions)**: MATLAB codes to compute IVIM-DKI parameter maps obtained from 
traditional IVIM-DKI and novel IVIM-DKI model with TV. 

    -'[hyModelTV.m](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/blob/main/IVIM_DKI_functions/hyModelTV.m)': Executes IVIM-DKI model with TV code using non-linear least square fitting with iterative total variation (TV) penalty function to perform spatial homogeneity on IVIM-DKI parameter reconstruction.
    
        -----------------------------------------------------------------------------------------------
        function[paraMap,resnorm,stats_roi]=hyModel(dwi,b,limit,initials,roi,stats,tvIter,alpha,const)
        -----------------------------------------------------------------------------------------------
        Input:
        dwi =     4D DWI data, MxNxSxB format where M and N are x and y, S is number of slices, 
                  and B is number of b-values 
        b =       b-values, must be a row matrix
        limit =   2x4 matrix with lower (1st row) and upper (2nd row) 
                  limits of all parameters in the order D, D*, f, and k
        initial = 1x4 matrix with initial values of all parameters 
                  in the order D, D*, f, and k
        roi =     Any region of interest, must be in 3D format MxNxS and logical
        stats =   If stats = true, then mean and std of ROI is calculated
        tvIter =  Define the number of TV iteration to be performed
        alpha =   Alpha is a small positive constant, range from 0 to 1
        const =   Const is a relaxation parameter, range from 0 to 1

        Output:
        paraMap =   IVIM-DKI parameters are saved as struct in the order D, D*, f, and k
        resnorm =   Voxelwise squared norm of the residual
        stats_roi = Mean and std of ROI are saved in [mean_roi, std_roi] format
                    in the order D, D*, f, and k

    -'[hymodel.m](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/blob/main/IVIM_DKI_functions/hyModel.m)': Executes traditional IVIM-DKI model code using ['allivimdki.m'](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/blob/main/IVIM_DKI_functions/allivimdki.m).
    
        -----------------------------------------------------------------------------
        function[paraMap,resnorm,stats_roi]=hyModel(dwi,b,limit,initials,roi,stats)
        -----------------------------------------------------------------------------
        Input:
        dwi =     4D DWI data, MxNxSxB format where M and N are x and y, S is number of slices, 
                  and B is number of b-values 
        b =       b-values, must be a row matrix
        limit =   2x4 matrix with lower (1st row) and upper (2nd row) 
                  limits of all parameters in the order D, D*, f, and k
        initial = 1x4 matrix with initial values of all parameters 
                  in the order D, D*, f, and k
        roi =     Any region of interest, must be in 3D format MxNxS and logical
        stats =   If stats = true, then mean and std of ROI is calculated

        Output:
        paraMap =   IVIM-DKI parameters are saved as struct in the order D, D*, f, and k
        resnorm =   Voxelwise squared norm of the residual
        stats_roi = Mean and std of ROI are saved in [mean_roi, std_roi] format
                    in the order D, D*, f, and k

  -'allivimdki.m': Function which contains IVIM-DKI model equation.
  
  -'[im2Y.m](https://www.mathworks.com/matlabcentral/fileexchange/65579-ivim-model-fitting)':Transforms functional image data (4D or 3D array) into data matrix VxM where V is the number of voxels and M is number of b-values.

  -'tv3d.m': Function which calculates 3D TV penalty. Please see ref. [1](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1002/mp.12520) for more details on the implementation of TV function.
  
**Note: Please execute 'main_program.m' in IVIM_DKI_functions folder.**

   - 'main_program.m': Example code which shows how to read or write the data and how to use the functions '[hyModelTV.m](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/blob/main/IVIM_DKI_functions/hyModelTV.m)' and '[hyModel.m](https://github.com/amitvmehndiratta/IVIM-DKI-MRMP2021/blob/main/IVIM_DKI_functions/hyModel.m)'.

If you use this function in research, please cite:

Malagi, A. V. et al. (2021). IVIM-DKI differentiation between prostate cancer and benign prostatic hyperplasia: comparison of 1.5T vs. 3T MRI. Magnetic Resonance Materials in Physics, Biology and Medicine.
doi:
