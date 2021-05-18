function f = monoexplog( Dm, b )  % Dm = ln(ADC)
% Monoexponential model
%   ADC parameter is unknown and free parameter
     C=exp(Dm);
     f=exp((-b).*C(1));
     f=f./f(1);
end

