%% AIC and AICc Calculations of all tissues for HY and HY+TV models
% AIC formula modified for least square fitting using residual sum of
% square

% Hybrid model
AIChymap_allb_tumor=AIChymap_allb.*tumor1;
AIChymap_allb_bph=AIChymap_allb.*bph1;
AIChymap_allb_pz=AIChymap_allb.*pz1;

% AICc

% Hybrid model
AICchymap_allb_tumor=AICchymap_allb.*tumor1;
AICchymap_allb_bph=AICchymap_allb.*bph1;
AICchymap_allb_pz=AICchymap_allb.*pz1;
%% Mean AIC and AICc values of all tissues
%AIC mean
mean_allAIChy=[nanmean( (AIChymap_allb_tumor(tumor1))) nanmean( (AIChymap_allb_bph(bph1))) nanmean( (AIChymap_allb_pz(pz1)))];
%nanstd
std_allAIChy=[nanstd( (AIChymap_allb_tumor(tumor1))) nanstd( (AIChymap_allb_bph(bph1))) nanstd( (AIChymap_allb_pz(pz1)))];
%AICc nanmean
mean_allAICchy=[nanmean( (AICchymap_allb_tumor(tumor1))) nanmean( (AICchymap_allb_bph(bph1))) nanmean( (AICchymap_allb_pz(pz1)))];
%nanstd
std_allAICchy=[nanstd( (AICchymap_allb_tumor(tumor1))) nanstd( (AICchymap_allb_bph(bph1))) nanstd( (AICchymap_allb_pz(pz1)))];
%% nanmean of AIC and AICc of HY&TV parameters
AIChytvmap_allb_tumor=AIChytvmap_allb.*tumor1;
AIChytvmap_allb_bph=AIChytvmap_allb.*bph1;
AIChytvmap_allb_pz=AIChytvmap_allb.*pz1;

AICchytvmap_allb_tumor=AICchytvmap_allb.*tumor1;
AICchytvmap_allb_bph=AICchytvmap_allb.*bph1;
AICchytvmap_allb_pz=AICchytvmap_allb.*pz1;

mean_allAIChytv=[nanmean( (AIChytvmap_allb_tumor(tumor1))) nanmean( (AIChytvmap_allb_bph(bph1))) nanmean( (AIChytvmap_allb_pz(pz1)))];
std_allAIChytv=[nanstd( (AIChytvmap_allb_tumor(tumor1))) nanstd( (AIChytvmap_allb_bph(bph1))) nanstd( (AIChytvmap_allb_pz(pz1)))];
mean_allAICchytv=[nanmean( (AICchytvmap_allb_tumor(tumor1))) nanmean( (AICchytvmap_allb_bph(bph1))) nanmean( (AICchytvmap_allb_pz(pz1)))];
std_allAICchytv=[nanstd( (AICchytvmap_allb_tumor(tumor1))) nanstd( (AICchytvmap_allb_bph(bph1))) nanstd( (AICchytvmap_allb_pz(pz1)))];
%% Saving all results in excel sheets
FeaturesRange1={'Tumor' 'BPH' 'Healthy PZ' 'Tumor' 'BPH' 'Healthy PZ'};
FeaturesRange2={'HY','HY+TV'};
FeaturesRangeValues1=[mean_allAIChy mean_allAIChytv];
FeaturesRangeValues2=[std_allAIChy std_allAIChytv];
FeaturesRangeValues3=[mean_allAICchy mean_allAICchytv];
FeaturesRangeValues4=[std_allAICchy std_allAICchytv];

xlswrite(strcat(pwd,'\aic_aicc.xlsx'),FeaturesRange1,'AIC','B1:G1');
xlswrite(strcat(pwd,'\aic_aicc.xlsx'),FeaturesRangeValues1,'AIC',strcat('B',num2str(n_mean),':G',num2str(n_mean)));
xlswrite(strcat(pwd,'\aic_aicc.xlsx'),FeaturesRange1,'AIC','I1:N1');
xlswrite(strcat(pwd,'\aic_aicc.xlsx'),FeaturesRangeValues2,'AIC',strcat('I',num2str(n_mean),':N',num2str(n_mean)));
xlswrite(strcat(pwd,'\aic_aicc.xlsx'),FeaturesRange1,'AICc','B1:G1');
xlswrite(strcat(pwd,'\aic_aicc.xlsx'),FeaturesRangeValues3,'AICc',strcat('B',num2str(n_mean),':G',num2str(n_mean)));
xlswrite(strcat(pwd,'\aic_aicc.xlsx'),FeaturesRange1,'AICc','I1:N1');
xlswrite(strcat(pwd,'\aic_aicc.xlsx'),FeaturesRangeValues4,'AICc',strcat('I',num2str(n_mean),':N',num2str(n_mean)));
