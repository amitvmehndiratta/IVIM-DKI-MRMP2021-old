function [s, gtv] = tv3d(img, flag)   
% Computes the total variation of the image.
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
    gtv = 0;
    [ux, uy, uz] = gradient(img); % ux <- x-gradient of the image
                                        % uy <- y-gradient of the image
    ds = sqrt(abs(ux).^2 + abs(uy).^2 + abs(uz).^2+ 1e-10); % Square of both ux and uy
    s = sum(ds(:));                     % Total Variation of the image; % V(x, y)                                    
                                        
    if(flag);
        uxx = ux./ds;                   % Normalized x-gradient
        uyy = uy./ds;                   % Normalized y-gradient
        uzz = uz./ds; 
        gtv = -divergence(uxx, uyy, uzz);    % Divergence x- and y- gradient
        gtv = gtv / norm( abs(gtv(:)), 'fro');
    end;
end
