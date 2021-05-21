%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IVIM-DKI differentiation between prostate cancer and benign %%%
%%%%%% prostatic hyperplasia: comparison of 1.5T vs. 3T MRI %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read patient data and ROI
% Select ivim file
[dwi_file,path]=uigetfile('*.nii.gz', 'Select IVIM-DKI data');
if isequal(dwi_file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(dwi_file)]);
end
dwiData=load_untouch_niigz(strcat(path,'\',dwi_file)); %load main image file
dwiData=double(dwiData.img);
dwiData=imrotate(dwiData,90);

% Select tumor file
[tumor_file,path]=uigetfile('*.nii.gz', 'Select tumor ROI data');
if isequal(tumor_file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(tumor_file)]);
end
tumor=load_untouch_niigz(strcat(path,'\',tumor_file)); %load main image file
tumor=double(tumor.img); 
tumor=logical(imrotate(tumor,90));

% Select bph file
[bph_file,path]=uigetfile('*.nii.gz', 'Select BPH ROI data');
if isequal(bph_file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(bph_file)]);
end
bph=load_untouch_niigz(strcat(path,'\',bph_file)); %load main image file
bph=double(bph.img);
bph=logical(imrotate(bph,90));

% Select pz file
[pz_file,path]=uigetfile('*.nii.gz', 'Select PZ ROI data');
if isequal(pz_file,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(pz_file)]);
end
pz=load_untouch_niigz(strcat(path,'\',pz_file)); %load main image file
pz=double(pz.img); 
pz=logical(imrotate(pz,90));

%% Define parameters 
%b-values 
 b=[0 25 50 75 100 150 200 500 800 1000 1250 1500 2000];
 
% To perform statistics
 stats=true;

% Select ROI
roi=tumor;

% Initial values, LB and UB of all pamameters
% Hybrid model and Hybrid model with TV
 D0_hy = [log(0.000729), log(0.013), log(0.1), log(1.19)]; 
 limit_hy= [log(0.0001), log(0.0001), log(0.001), log(0.01);...
    log(0.05), log(0.5), log(1), log(3)];

%TV parameters(alpha, const and TV iteration)
 alpha=0.005;
 const=0.99;
 TViter=10;

%% Codes to run
 [paraMapHY,resnorm,stats_roiHY] = hyModel(dwiData,b,limit_hy,D0_hy,tumor,stats);
 
 [paraMapHYTV,resnormTV,stats_roiHYTV] = hyModelTV(dwiData,b,limit_hy,D0_hy,tumor,...
    stats,TViter,alpha,const);

%% Save parameter maps in compressed nifti format
mkdir(strcat(path,'Output'))

% Save parameters from HY model
save_untouch_niigz(imrotate(paraMap.DmapHY,-90),...
    strcat(path,'Output\Dmap_HY.nii.gz'),strcat(path,dwi_file))
save_untouch_niigz(imrotate(paraMap.DpmapHY,-90),...
    strcat(path,'Output\Dpmap_HY.nii.gz'),strcat(path,dwi_file))
save_untouch_niigz(imrotate(paraMap.fmapHY,-90),...
    strcat(path,'Output\fmap_HY.nii.gz'),strcat(path,dwi_file))
save_untouch_niigz(imrotate(paraMap.kmapHY,-90),...
    strcat(path,'Output\kmap_HY.nii.gz'),strcat(path,dwi_file))

% Save parameters from HY model with TV
save_untouch_niigz(imrotate(paraMapTV.DmapHYtv,-90),...
    strcat(path,'Output\Dmap_HYtv.nii.gz'),strcat(path,dwi_file))
save_untouch_niigz(imrotate(paraMapTV.DpmapHYtv,-90),...
    strcat(path,'Output\Dpmap_HYtv.nii.gz'),strcat(path,dwi_file))
save_untouch_niigz(imrotate(paraMapTV.fmapHYtv,-90),...
    strcat(path,'Output\fmap_HYtv.nii.gz'),strcat(path,dwi_file))
save_untouch_niigz(imrotate(paraMapTV.kmapHYtv,-90),...
    strcat(path,'Output\kmap_HYtv.nii.gz'),strcat(path,dwi_file))
