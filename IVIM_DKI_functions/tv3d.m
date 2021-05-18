function [s, gtv] = tv3d(img, flag)   
%% Computes the total variation of the image.
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
