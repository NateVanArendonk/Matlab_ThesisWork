function [cp, cg] = getcp(f,h)
% [cp, cg] = getc(f,h)
g = 9.81;
k = getk(f,h);
cp = sqrt(g./k.*tanh(k.*h));
%cg = 0.5*(g*tanh(k*h)+g*(k*h)*(sech(k*h)^2))/sqrt(g*k*tanh(k*h));
% cg = cp./2*(1+2*k.*h/sinh(2*k.*h)); %Kundu & Cohen p240
cg = cp./(2*(1+2*k.*h./sinh(2*k.*h)));

end

function k = getk(f,h) 
%k = getk(f,h)
%Credit F. Fedderson
% returns the wavenumber of the gravity wave
% dispersion relation, by using newtons method
% the initial guess will be the shallow water wavenumber
ogf = f;
ogh = h;
omega = 2*pi*f;
g = 9.81;
k = omega./sqrt(g*h);
f = g*k.*tanh(k.*h) - omega.^2;
count = 0;
while max(abs(f))>1e-6
    dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
    k = k - f./dfdk;
    f = g*k.*tanh(k.*h) - omega.^2;
%     if count > 1000
%         fprintf('It is Hung\n')
%         fprintf('F: %.2f and h: %.2f',ogf,ogh);
%     end
%     count = count+1;
end
end