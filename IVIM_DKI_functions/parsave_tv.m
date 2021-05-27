function parsave_tv(fname, x,y,z,k,res)
% Saving variables from parallel loop 'parfor'
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
save(fname, 'x','y','z','k','res')
end