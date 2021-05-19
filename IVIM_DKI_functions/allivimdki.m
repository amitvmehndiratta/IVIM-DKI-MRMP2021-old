function ydata = allivimdki( Dp, xdata )   
% IVIM-DKI hybrid model 
% Four free parameter, D, Dp, f, and k
   para=exp(Dp);
   ydata=para(3)*exp( -xdata*(para(2)) ) + (1-para(3))...
       *exp(-xdata*para(1)+((((xdata*para(1))).^2)*para(4))/6 );
end