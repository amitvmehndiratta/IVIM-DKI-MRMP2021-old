%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IVIM-DKI differentiation between prostate cancer and benign %%%
%%%%%% prostatic hyperplasia: comparison of 1.5T vs. 3T MRI %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read patient data
% Select ivim file
[dwi_file,path]=uigetfile('ivim13b.nii.gz');
if isequal(dwi_file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(dwi_file)]);
end
allslice=load_untouch_niigz(strcat(path,'\',dwi_file)); %load main image file
allslice=double(allslice.img); %double the image for further easy processing
allslice3=imrotate(allslice,90);

% Select tumor file
[tumor_file,path]=uigetfile('tumor.nii.gz');
if isequal(tumor_file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(tumor_file)]);
end
tumor=load_untouch_niigz(strcat(path,'\',tumor_file)); %load main image file
tumor=double(tumor.img); %double the image for further easy processing
tumor1=logical(imrotate(tumor,90));

% Select bph file
[bph_file,path]=uigetfile('bph.nii.gz');
if isequal(bph_file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(bph_file)]);
end
bph=load_untouch_niigz(strcat(path,'\',bph_file)); %load main image file
bph=double(bph.img); %double the image for further easy processing
bph1=logical(imrotate(bph,90));

% Select pz file
[pz_file,path]=uigetfile('pz.nii.gz');
if isequal(pz_file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(pz_file)]);
end
pz=load_untouch_niigz(strcat(path,'\',pz_file)); %load main image file
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
%% Save workspace
save (strcat(path,'allmodel_results.mat'))