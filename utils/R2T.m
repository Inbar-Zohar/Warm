function T = R2T(R)
% TERMISTOR NAME:
% Resistance @ room temperature: 10.5~12.5 kOhm

% % old:
% p = [  250.7, 5551   , 7240   ,-3672   , -693.4 ];
% q = [    1  ,   17.71,   17.33,  -10.65,   -1.31];
% T = polyval(p,R)./polyval(q,R);

p = [3.661e5,3.035e8,1.472e10];
q = [1.000  ,4954.0 ,1.824e6 ,5.145e7];
T = polyval(p,R)./polyval(q,R);

end