%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IVIM-DKI differentiation between prostate cancer and benign %%%
%%%%%% prostatic hyperplasia: comparison of 1.5T vs. 3T MRI %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read patient data
patient='p1';
mri='1_5t';
%Read patient data
path = pwd;% path of working directory
dwi_filename = 'ivim13b.nii.gz'; %image file names
tumor_filename = 'tumor.nii.gz';
bph_filename = 'bph.nii.gz';
pz_filename = 'pz.nii.gz';

allslice=load_untouch_niigz(strcat(path,'\',dwi_filename)); %load main image file
allslice=double(allslice.img); %double the image for further easy processing
allslice3=imrotate(allslice,90);

tumor=load_untouch_niigz(strcat(path,'\',tumor_filename)); %load main image file
tumor=double(tumor.img); %double the image for further easy processing
tumor1=logical(imrotate(tumor,90));

bph=load_untouch_niigz(strcat(path,'\',bph_filename)); %load main image file
bph=double(bph.img); %double the image for further easy processing
bph1=logical(imrotate(bph,90));

pz=load_untouch_niigz(strcat(path,'\',pz_filename)); %load main image file
pz=double(pz.img); %double the image for further easy processing
pz1=logical(imrotate(pz,90));

%% Constants 
%b-values 
 b=[0 25 50 75 100 150 200 500 800 1000 1250 1500 2000];
 b_mono=[0 500 800 1000];
%TV parameters(alpha, beta and TV iteration)
alpha1=0.005;
const=0.99;
iter=10;
%For AIC
 num=length(b);
 k_hy=4;
 n_mean=2;
%% Codes to run
hybridmodel
hybrid_TVmodel
aic_aicccal
%% Save all results as mat file
save (strcat('allmodel_results',mri,patient,'.mat'))