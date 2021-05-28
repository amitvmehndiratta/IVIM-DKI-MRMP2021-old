function [mean_roi, std_roi] = stats_roi(paraMap,roi)
%-----------------------------------------------------------------------------
% function[mean_roi, std_roi] = stats_roi(parameterMap,roi,stats)
%-----------------------------------------------------------------------------
% Description: Calculates mean and standard deviation of ROI provided 
%              for IVIM-DKI parameters
%
% Input:
% paraMap =   IVIM-DKI parameters saved as struct in the order D, D*, f, and k (output generated from hyModelTV.m)
% roi =       ROI mask in 2D or 3D logical matrix format

% Output:
% mean_roi = Average of IVIM-DKI parameters for ROI specified and 
%            are saved as struct in the order D, D*, f, and k
% std_roi = Standard deviation of IVIM-DKI parameters for ROI specified and 
%            are saved as struct in the order D, D*, f, and k
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
%% Statistics of ROI provided
    fNames=fieldnames(paraMap);
    for i=1:length(fNames)
        mean_roi.(string(fNames(i)))=mean(paraMap.(string(fNames(i)))(roi));
        std_roi.(string(fNames(i)))=std(paraMap.(string(fNames(i)))(roi));
    end

end
