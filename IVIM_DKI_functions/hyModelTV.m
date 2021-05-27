%% IVIM-DKI Model fitting with TV penalty function

function [paraMap,resnorm] = hyModelTV(dwi,b,ilimit,initial,tvIter,alpha,beta)
%-----------------------------------------------------------------------------------------------
% function[paraMap,resnorm,stats_roi]=hyModel(dwi,b,ilimit,initials,roi,stats,tvIter,alpha,beta)
%-----------------------------------------------------------------------------------------------
% Description: Executes IVIM-DKI model with TV code using non-linear least 
%              square optimization with iterative total variation (TV) 
%              penalty function to evaluate IVIM-DKI parameters.
% Input:
% dwi =     4D DWI data, MxNxSxB format where M and N are x and y, S is number of slices, 
%           and B is number of b-values 
% b =       b-values, must be a row matrix
% ilimit =   2x4 matrix with lower (1st row) and upper (2nd row) 
%           ilimits of all parameters in the order D, D*, f, and k
% initial = 1x4 matrix with initial values of all parameters 
%           in the order D, D*, f, and k
% tvIter =  Define the number of TV iteration to be performed
% alpha =   Alpha is a small positive betaant, range from 0 to 1
% beta =   beta is a relaxation parameter, range from 0 to 1

% Output:
% paraMap =   IVIM-DKI parameters are saved as struct in the order D, D*, f, and k
% resnorm =   Voxelwise squared norm of the residual
%
% Copyright © 2021 IIT Delhi, India.
% Disclaimer: This project can be used only for research purposes. Authors 
% are not liable for any clinical use of it, authors could not be held responsible.
% For queries or suggestions about this package, kindly contact Dr. Amit Mehndiratta, 
% IIT Delhi, India. E-mail: amehndiratta@cbme.iitd.ac.in, amit.mehndiratta@keble.oxon.org.
%
% If you use this function in research, please cite:
%
% Malagi, A. V. et al. (2021). IVIM-DKI differentiation between prostate cancer 
% and benign prostatic hyperplasia: comparison of 1.5T vs. 3T MRI. 
% Magnetic Resonance Materials in Physics, Biology and Medicine.
% doi:10.1007/s10334-021-00932-1    
%% Check the inputs given
if length(size(dwi))~=4
    error('DWI image must be a 4D matrix');
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
 %% Initialization of parameter maps
DmapHYtv_allb=zeros(row,col,totalslice);
DpmapHYtv_allb=zeros(row,col,totalslice);
fmapHYtv_allb=zeros(row,col,totalslice);
kmapHYtv_allb=zeros(row,col,totalslice);
resnorm=zeros(row,col,totalslice);
%% IVIM-DKI with TV fitting
start=tic;
%% IVIM-DKI with TV and parallel computing
% lsqcurvefitting & TV reduction for hy fitting
parfor smoothn=1:tvIter % TV iteration
    % Vectorization
    ydata_ivimdki=im2Y(dwiSignal,mask);
    [vox,~]=size(ydata_ivimdki);
    % Initialization of parameter maps
    Dbi_tv= zeros(row,col,totalslice);
    Dpbi_tv= zeros(row,col,totalslice);
    fbi_tv= zeros(row,col,totalslice);
    kbi_tv= zeros(row,col,totalslice);
    
    Dbi_tv(mask)=initial(1);
    Dpbi_tv(mask)=initial(2);
    fbi_tv(mask)=initial(3);
    kbi_tv(mask)=initial(4);

    Dmaptv1=zeros(vox,1);
    Dpmaptv1=zeros(vox,1);
    fmaptv1=zeros(vox,1);
    kmaptv1=zeros(vox,1);
    resHYtv_allb1=zeros(vox,1);

    Dmaptv=zeros(row,col,totalslice);
    Dpmaptv=zeros(row,col,totalslice);
    fmaptv=zeros(row,col,totalslice);
    kmaptv=zeros(row,col,totalslice);
     
    fprintf(' Processing TV iter: %d \n',smoothn);
    Dtv=Dbi_tv(mask);
    Dptv=Dpbi_tv(mask);
    ftv=fbi_tv(mask);
    ktv=kbi_tv(mask);

    for slice=1:vox
        options = optimset('MaxIter',1, 'Display','off');  %  lsqcurvefitting for 1 iteration
        ydata_13b=squeeze(ydata_ivimdki(slice,:));
        D0_hy2 = [Dtv(slice), Dptv(slice), ftv(slice), ktv(slice)]; 
        [HYtvDDpfk, resnormHYtv] = lsqcurvefit(@allivimdki,D0_hy2,b,ydata_13b,ilimit(1,:),ilimit(2,:),options);
        % assign D, D*, f and K value per voxel
        Dmaptv1(slice,:) = HYtvDDpfk(1);
        Dpmaptv1(slice,:) = HYtvDDpfk(2);
        fmaptv1(slice,:) = HYtvDDpfk(3);
        kmaptv1(slice,:) = HYtvDDpfk(4);
        resHYtv_allb1(slice,:) = resnormHYtv;
    end

    Dmaptv(mask )=Dmaptv1;
    Dpmaptv(mask )=Dpmaptv1;
    fmaptv(mask )=fmaptv1;
    kmaptv(mask )=kmaptv1;

    % TV reduction for D diffusion coeff
    D_tv = Dmaptv;
    d1b = norm(Dbi_tv(:) - D_tv(:) + 1e-10, 'fro');
    for nDtv= 1:2000
        [s1,gtv1] = tv3d(Dmaptv,1);
        Dmaptv = Dmaptv - alpha*d1b*gtv1;
        d2b = norm(Dmaptv(:) - D_tv(:) + 1e-10, 'fro');
        if(d2b > (beta*d1b))
            break;
        end
    end   % TV reduction ends

    Dbi_tv = Dmaptv;

% TV reduction for Dp perfusion coeff

    Dp_tv = Dpmaptv;
    d1b = norm(Dpbi_tv(:) - Dp_tv(:) + 1e-10, 'fro');
    for nDptv= 1:2000
        [s2,gtv2] = tv3d(Dpmaptv,1);
        Dpmaptv = Dpmaptv - alpha*d1b*gtv2;
        d2b = norm(Dpmaptv(:) - Dp_tv(:) + 1e-10, 'fro');
        if(d2b > (beta*d1b))
           break;
        end
    end   % TV reduction ends

    Dpbi_tv = Dpmaptv;


    % TV reduction for f perfusion fraction

    f_tv = fmaptv;
    d1c = norm(fbi_tv(:) - f_tv(:) + 1e-10, 'fro');

    for nftv= 1:2000
        [s3,gtv3] = tv3d(fmaptv,1);
        fmaptv = fmaptv - alpha*d1c*gtv3;
        d2c = norm(fmaptv(:) - f_tv(:) + 1e-10, 'fro');
        if(d2c > (beta*d1c))
            break;
        end
     end   % TV reduction ends

    fbi_tv = fmaptv;

    % TV reduction for k kurtosis

    k_tv = kmaptv;
    d1d = norm(kbi_tv(:) - k_tv(:) + 1e-10, 'fro');

    for nktv= 1:2000
        [s4,gtv4] = tv3d(kmaptv,1);
        kmaptv = kmaptv - alpha*d1d*gtv4;
        d2d = norm(kmaptv(:) - k_tv(:) + 1e-10, 'fro');
        if(d2d > (beta*d1d))
            break;
        end
    end   % TV reduction ends

    kbi_tv = kmaptv;
    parsave_tv(sprintf('output_tv%d.mat', smoothn), D_tv, Dp_tv, f_tv,k_tv,resHYtv_allb1);
end
totalTime_TV=toc(start);
 fprintf(strcat('Total time taken by novel HY model with TV:',...
    num2str(round(floor(totalTime_TV/60),0)),'.',num2str(round(rem(totalTime_TV,60),0)),' minutes\n'))
%% Saving Parameter maps
load(strcat(pwd,'\output_tv10.mat')) % Take the last iteration output from parloop
DmapHYtv_allb(mask )=exp(x(mask ));
DpmapHYtv_allb(mask )=exp(y(mask ));
fmapHYtv_allb(mask )=exp(z(mask ));
kmapHYtv_allb(mask )=exp(k(mask ));
resnorm(mask)=exp(res);
paraMap.DmapHYtv=DmapHYtv_allb;
paraMap.DpmapHYtv=DpmapHYtv_allb;
paraMap.fmapHYtv=fmapHYtv_allb;
paraMap.kmapHYtv=kmapHYtv_allb;
% Delete temp files from directory
for i=1:tvIter
    delete(strcat(pwd,'\output_tv',num2str(i),'.mat'))
end

end