%% IVIM-DKI analysis using Hybrid model
 [~,~,slice,~] = size(allslice3(:,:,:,:));
% Body mask generation (eliminating background noise)
for s=1:slice
        thresh =10;
        thrI=imbinarize(allslice3(:,:,s,1),thresh);
        se = strel('disk',3);
        x=imclose(thrI,se);
        se1 = strel('disk',35);
        x1=imerode(x,se1);
        se2 = strel('disk',19);
        x2=imdilate(x1,se2);
        allslice1(:,:,s)=allslice3(:,:,s,1).*x2;
end
mask=logical(allslice1);
allslice11=allslice3.*mask;
%% Initial values, LB and UB of all pamameters
%Monoexp model
D0_mono = log(0.000729); 
lb_mono= log(0.0001);
ub_mono = log(0.05);
% Hybrid model
D0_hy = [log(0.000729), log(0.013), log(0.1), log(1.19)]; 
lb_ivimdki= [log(0.0001), log(0.0001), log(0.001), log(0.01)];
ub_ivimdki = [log(0.05), log(0.5), log(1), log(3)];
 %[row,col,totalslice,totalb] = size(allslice11(:,:,:,:));
 %% Normalization of IVIM-DKI data
Smain1=allslice11./allslice11(:,:,:,1);
Smain1(isnan(Smain1))=0;
Smain1(isinf(Smain1))=0;
Smain1_mono=Smain1(:,:,:,[1,8:10]);
[row1,col,totalslice,~] = size(Smain1(:,:,:,:));
ydata_ivimdki=im2Y(Smain1,mask); % Vectorization
ydata_mono=im2Y(Smain1_mono,mask);
[row,totalb]=size(ydata_ivimdki);
%% Initialization of maps
%Parameter maps
ADCmap_allb=zeros(row1,col,totalslice);
Dmaphy_allb=zeros(row1,col,totalslice);
Dpmaphy_allb=zeros(row1,col,totalslice);
fmaphy_allb=zeros(row1,col,totalslice);
kmaphy_allb=zeros(row1,col,totalslice);
resmaphy_allb=zeros(row1,col,totalslice);
AIChymap_allb=zeros(row1,col,totalslice);
AICchymap_allb=zeros(row1,col,totalslice);
    %Temp parameter for parallel loop
ADC_allb=zeros(row,1);
dhy_allb=zeros(row,1);
dphy_allb=zeros(row,1);
fhy_allb=zeros(row,1);
khy_allb=zeros(row,1);
reshy_allb=zeros(row,1);
AIChy_allb=zeros(row,1);
AICchy_allb= zeros(row,1);
%% Model fitting
% profile on -history
tic
for slice= 1:row %parallel computing
fprintf(' Processing voxel: %d \n',slice);   
options = optimset('MaxIter',[], 'Display','off');   
       ydata13b_mono=squeeze(ydata_mono(slice,:));
       ydata13b=squeeze(ydata_ivimdki(slice,:));
       %ydata_biD1013b=ydata_biD1013b';
%% Monoexponential model fit using all b-values
       [monoADC, resnormmono, res_mono] = lsqcurvefit(@monoexplog,D0_mono,b_mono,ydata13b_mono,lb_mono,ub_mono,options);
% assign ADC value per voxel
       ADC_allb(slice,:) = monoADC;

%% Hybrid model fit using all b-values
       [hyDDpfk, resnormhy, res_hy] = lsqcurvefit(@allivimdki,D0_hy,b,ydata13b,lb_ivimdki,ub_ivimdki,options);
% assign D, D*, f and k value per voxel
       dhy_allb(slice,:) = hyDDpfk(1);
       dphy_allb(slice,:) = hyDDpfk(2);
       fhy_allb(slice,:) = hyDDpfk(3);
       khy_allb(slice,:) = hyDDpfk(4);
       reshy_allb(slice,:)=resnormhy;
       [AIChy_allb(slice,:),AICchy_allb(slice,:)]=AIC(exp(resnormhy),length(res_hy),k_hy);
       
end

%% Saving Parameter maps
ADCmap_allb(mask)=exp(ADC_allb);
Dmaphy_allb(mask)=exp(dhy_allb);
Dpmaphy_allb(mask)=exp(dphy_allb);
fmaphy_allb(mask)=exp(fhy_allb);
kmaphy_allb(mask)=exp(khy_allb);
resmaphy_allb(mask)=exp(reshy_allb);
AIChymap_allb(mask)=(AIChy_allb);
AICchymap_allb(mask)=(AICchy_allb);

total_time_allmod=toc;