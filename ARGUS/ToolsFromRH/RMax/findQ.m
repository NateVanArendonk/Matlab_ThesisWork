function Q = findQ(params, Hs, dt)
%   Q = findQ(params, Hs, dt)
%
% estimate a wave-height dependent process error, Q using specified params.
% Selected form is Q = params(1) + params(2) * (max(Hs-1,0))^2.  i.e. it is
% dependent of wave height but floors for Hs < 1.  Final units must be
% variance (m^2) so params must yield terms in m^2/deltaT, where deltaT is
% in days.  

Q = params(1) + params(2) * (max(Hs-1,0))^2;
Q = Q * dt;