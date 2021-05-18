%% IVIM-DKI Model fitting with TV penalty function

% Initialization for TV pamameter fitting
lb_ivimdkitv= [log(0.0001), log(0.0001), log(0.001), log(0.01)];
ub_ivimdkitv = [log(0.05), log(0.5), log(1), log(3)];
 %[row,col,totalslice,totalb] = size(allslice11(:,:,:,:));

[row1,col,totalslice,~] = size(Smain1(:,:,:,:));
D_initial= log(0.000729);
Dp_initial = log(0.013);
f_initial = log(0.1);
k_initial = log(1.19);
% Initialization of parameters
 %% Initialization of parameter maps
Dmaphytv_allb=zeros(row1,col,totalslice);
Dpmaphytv_allb=zeros(row1,col,totalslice);
fmaphytv_allb=zeros(row1,col,totalslice);
kmaphytv_allb=zeros(row1,col,totalslice);
reshytvmap_allb=zeros(row1,col,totalslice);
AIChytvmap_allb=zeros(row1,col,totalslice);
AICchytvmap_allb=zeros(row1,col,totalslice);
%% IVIM-DKI with TV fitting
% profile on -history

tic
%% IVIM-DKI with TV and parallel computing
% lsqcurvefitting & TV reduction for hy fitting
parfor smoothn=1:iter % TV iteration
    ydata_ivimdki=im2Y(Smain1,mask);
    [row,totalb]=size(ydata_ivimdki);
Dbi_tv= zeros(row1,col,totalslice);
Dpbi_tv= zeros(row1,col,totalslice);
fbi_tv= zeros(row1,col,totalslice);
kbi_tv= zeros(row1,col,totalslice);
    
Dbi_tv(mask)=D_initial;
Dpbi_tv(mask)=Dp_initial;
fbi_tv(mask)=f_initial;
kbi_tv(mask)=k_initial;

Dmaptv1=zeros(row,1);
Dpmaptv1=zeros(row,1);
fmaptv1=zeros(row,1);
kmaptv1=zeros(row,1);
reshytv_allb1=zeros(row,1);
AIChytv_allb1=zeros(row,1);
AICchytv_allb1= zeros(row,1);

Dmaptv=zeros(row1,col,totalslice);
Dpmaptv=zeros(row1,col,totalslice);
fmaptv=zeros(row1,col,totalslice);
kmaptv=zeros(row1,col,totalslice);
     
    fprintf(' Processing Tv iter: %d \n',smoothn);
    Dtv=Dbi_tv(mask);
    Dptv=Dpbi_tv(mask);
    ftv=fbi_tv(mask);
    ktv=kbi_tv(mask);

for slice=1:row
options = optimset('MaxIter',1, 'Display','off');  %  lsqcurvefitting for 1 iteration
ydata_biD1013b=squeeze(ydata_ivimdki(slice,:));
       D0_hy2 = [Dtv(slice), Dptv(slice), ftv(slice), ktv(slice)]; 
       [hytvDDpfk, resnormhytv, res_hytv] = lsqcurvefit(@allivimdki,D0_hy2,b,ydata_biD1013b,lb_ivimdkitv,ub_ivimdkitv,options);
  
   % assign D, D*, f and K value per voxel
      Dmaptv1(slice,:) = hytvDDpfk(1);
      Dpmaptv1(slice,:) = hytvDDpfk(2);
      fmaptv1(slice,:) = hytvDDpfk(3);
      kmaptv1(slice,:) = hytvDDpfk(4);
      reshytv_allb1(slice,:) = resnormhytv;
      [AIChytv_allb1(slice,:),AICchytv_allb1(slice,:)]=AIC(exp(resnormhytv),length(res_hytv),k_hy);
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
    Dmaptv = Dmaptv - alpha1*d1b*gtv1;
    d2b = norm(Dmaptv(:) - D_tv(:) + 1e-10, 'fro');
     if(d2b > (const*d1b))
        break;
     end
  end   % TV reduction ends

Dbi_tv = Dmaptv;

% TV reduction for Dp perfusion coeff

Dp_tv = Dpmaptv;
d1b = norm(Dpbi_tv(:) - Dp_tv(:) + 1e-10, 'fro');
  for nDptv= 1:2000
   
    [s2,gtv2] = tv3d(Dpmaptv,1);
    Dpmaptv = Dpmaptv - alpha1*d1b*gtv2;
    d2b = norm(Dpmaptv(:) - Dp_tv(:) + 1e-10, 'fro');
     if(d2b > (const*d1b))
        break;
     end
  end   % TV reduction ends

Dpbi_tv = Dpmaptv;


% TV reduction for f perfusion fraction

f_tv = fmaptv;
d1c = norm(fbi_tv(:) - f_tv(:) + 1e-10, 'fro');

  for nftv= 1:2000
   
    [s3,gtv3] = tv3d(fmaptv,1);
    fmaptv = fmaptv - alpha1*d1c*gtv3;
    d2c = norm(fmaptv(:) - f_tv(:) + 1e-10, 'fro');
     if(d2c > (const*d1c))
        break;
     end
  end   % TV reduction ends

fbi_tv = fmaptv;

% TV reduction for k kurtosis

k_tv = kmaptv;
d1d = norm(kbi_tv(:) - k_tv(:) + 1e-10, 'fro');

  for nktv= 1:2000
   
    [s4,gtv4] = tv3d(kmaptv,1);
    kmaptv = kmaptv - alpha1*d1d*gtv4;
    d2d = norm(kmaptv(:) - k_tv(:) + 1e-10, 'fro');
     if(d2d > (const*d1d))
        break;
     end
  end   % TV reduction ends

kbi_tv = kmaptv;
parsave_tv(sprintf('output_tv%d.mat', smoothn), D_tv, Dp_tv, f_tv,k_tv,reshytv_allb1,AIChytv_allb1,AICchytv_allb1);
end
total_time_TV=toc;
load('output_tv10.mat') % Take the last iteration output
Dmaphytv_allb(mask )=exp(x(mask ));
Dpmaphytv_allb(mask )=exp(y(mask ));
fmaphytv_allb(mask )=exp(z(mask ));
kmaphytv_allb(mask )=exp(k(mask ));
reshytvmap_allb=exp(ar);
AIChytvmap_allb(mask )=(br);
AICchytvmap_allb(mask )=(cr);