function R2 = StockdonEquation(H0, fp, beta)
%   R2 = StockdonEquation(H0, fp, beta);
%
% This function estimates the 2% runup exceedance, R2, as a function of the
% deep water wave height, H0 (m), the deep water frequency, fp (Hz) and the
% foreshore beach slope, beta.  The equations are taken from Stockdon and
% Holman, 2006.

g = 9.8;
L0 = g/(2*pi*fp*fp);
Ir = beta/(sqrt(H0/L0));
HL = H0*L0;

R2 = repmat(0.043*sqrt(HL), size(beta));    % default to dissipative
if Ir>=0.3
    R2 = 1.1*(0.35*beta*sqrt(HL) + sqrt(HL*(0.563*beta.*beta+0.004))/2);
end
