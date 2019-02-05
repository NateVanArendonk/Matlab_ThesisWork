function k = getk(f,h) 
%k = getk(f,h)
%Credit F. Fedderson
% returns the wavenumber of the gravity wave
% dispersion relation, by using newtons method
% the initial guess will be the shallow water wavenumber
omega = 2*pi*f;
g = 9.81;
k = omega./sqrt(g*h);
f = g*k.*tanh(k.*h) - omega.^2;
while max(abs(f))>1e-6
    dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
    k = k - f./dfdk;
    f = g*k.*tanh(k.*h) - omega.^2;
end
end