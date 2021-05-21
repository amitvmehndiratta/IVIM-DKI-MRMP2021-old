%% IVIM-DKI analysis using Hybrid model
% Simultaneously generates IVIM-DKI parameters
function [paraMap,resnorm,stats_roi] = hyModel(dwi,b,limit,initials,roi,stats)
%-----------------------------------------------------------------------------
% function[paraMap,resnorm,stats_roi]=hyModel(dwi,b,limit,initials,roi,stats)
%-----------------------------------------------------------------------------
% Input:
% dwi =     4D DWI data, MxNxSxB format where M and N are x and y, S is number of slices, 
%           and B is number of b-values 
% b =       b-values, must be a row matrix
% limit =   2x4 matrix with lower (1st row) and upper (2nd row) 
%           limits of all parameters in the order D, D*, f, and k
% initial = 1x4 matrix with initial values of all parameters 
%           in the order D, D*, f, and k
% roi =     Any region of interest, must be in 3D format MxNxS and logical
% stats =   If stats = true, then mean and std of ROI is calculated

% Output:
% paraMap =   IVIM-DKI parameters are saved as struct in the order D, D*, f, and k
% resnorm =   Voxelwise squared norm of the residual
% stats_roi = Mean and std of ROI are saved in [mean_roi, std_roi] format
%             in the order D, D*, f, and k
%
% If you use this function in research, please cite:
%
% Malagi, A. V. et al. (2021). IVIM-DKI differentiation between prostate cancer 
% and benign prostatic hyperplasia: comparison of 1.5T vs. 3T MRI. 
% Magnetic Resonance Materials in Physics, Biology and Medicine.
% doi:
%% Check the inputs given
if length(size(dwi))~=4
    error('DWI image must be a 4D matrix');
end
if length(size(roi))~=3
    error('ROI data must be a 3D matrix');
end
if ~islogical(roi)
    error('ROI data must be logical');
end
if ~isrow(b)
    b = b';
end
%% Body mask generation (eliminating background noise)
for s=1:size(dwi,3)
        thresh =10;
        thrI=imbinarize(dwi(:,:,s,1),thresh);
        se = strel('disk',3);
        x=imclose(thrI,se);
        se1 = strel('disk',35);
        x1=imerode(x,se1);
        se2 = strel('disk',19);
        x2=imdilate(x1,se2);
        allslice1(:,:,s)=dwi(:,:,s,1).*x2;
end
mask=logical(allslice1);
dwiMasked=dwi.*mask;
%% Normalization of IVIM-DKI data
dwiSignal=dwiMasked./dwiMasked(:,:,:,1);
dwiSignal(isnan(dwiSignal))=0;
dwiSignal(isinf(dwiSignal))=0;
[row,col,totalslice,~] = size(dwiSignal);
%% Vectorization
ydata_ivimdki=im2Y(dwiSignal,mask); 
[vox,~]=size(ydata_ivimdki);
%% Initialization of maps
%Parameter maps
ADCmap_allb=zeros(row,col,totalslice);
DmapHY_allb=zeros(row,col,totalslice);
DpmapHY_allb=zeros(row,col,totalslice);
fmapHY_allb=zeros(row,col,totalslice);
kmapHY_allb=zeros(row,col,totalslice);
resnorm=zeros(row,col,totalslice);
%Temp parameters
dHY_allb=zeros(vox,1);
dpHY_allb=zeros(vox,1);
fHY_allb=zeros(vox,1);
kHY_allb=zeros(vox,1);
resHY_allb=zeros(vox,1);
%% Model fitting
fprintf('Processing IVIM-DKI data...\n')
start=tic;
for slice= 1:vox   
       options = optimset('MaxIter',[], 'Display','off');   
       ydata13b=squeeze(ydata_ivimdki(slice,:));
%% Hybrid model fit using all b-values
       [HYDDpfk, resnormHY] = lsqcurvefit(@allivimdki,initials,b,ydata13b,limit(1,:),limit(2,:),options);
% Assign D, D*, f and k value per voxel
       dHY_allb(slice,:) = HYDDpfk(1);
       dpHY_allb(slice,:) = HYDDpfk(2);
       fHY_allb(slice,:) = HYDDpfk(3);
       kHY_allb(slice,:) = HYDDpfk(4);
       resHY_allb(slice,:)=resnormHY;
end
totalTime_hy=toc(start);
 fprintf(strcat('Total time taken by traditional HY model:',...
    num2str(round(floor(totalTime_hy/60),0)),'.',num2str(round(rem(totalTime_hy,60),0)),' minutes\n'))
%% Saving Parameter maps
DmapHY_allb(mask)=exp(dHY_allb);
DpmapHY_allb(mask)=exp(dpHY_allb);
fmapHY_allb(mask)=exp(fHY_allb);
kmapHY_allb(mask)=exp(kHY_allb);
resnorm(mask)=exp(resHY_allb);
paraMap.DmapHY=DmapHY_allb;
paraMap.DpmapHY=DpmapHY_allb;
paraMap.fmapHY=fmapHY_allb;
paraMap.kmapHY=kmapHY_allb;
%% Statistics of ROI provided
if stats==1
    % Calculating mean of ROI given
    meanROI_Dmap=mean(DmapHY_allb(roi));
    meanROI_Dpmap=mean(DpmapHY_allb(roi));
    meanROI_fmap=mean(fmapHY_allb(roi));
    meanROI_kmap=mean(kmapHY_allb(roi));
    % Calculating standard deviation of ROI given
    stdROI_Dmap=std(DmapHY_allb(roi));
    stdROI_Dpmap=std(DpmapHY_allb(roi));
    stdROI_fmap=std(fmapHY_allb(roi));
    stdROI_kmap=std(kmapHY_allb(roi));
    %Output saving mean and std of ROI
    stats_roi.D=[meanROI_Dmap,stdROI_Dmap];
    stats_roi.Dp=[meanROI_Dpmap,stdROI_Dpmap];
    stats_roi.f=[meanROI_fmap,stdROI_fmap];
    stats_roi.k=[meanROI_kmap,stdROI_kmap]; 
end
end
