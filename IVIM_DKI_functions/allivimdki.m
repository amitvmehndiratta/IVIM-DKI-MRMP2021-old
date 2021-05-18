function v= allivimdki( Dp, xdata )   
% IVIM-DKI hybrid model 
% D, D*, f and K are calculated by non-linear curve fitting
% Four free parameter
   Bk=exp(Dp);
   v=Bk(3)*exp( -xdata*(Bk(2)) ) + (1-Bk(3))*exp( -xdata*Bk(1)+((((xdata*Bk(1))).^2)*Bk(4))/6 );
end