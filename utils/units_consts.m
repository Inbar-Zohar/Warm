function [U,P,C] = units_consts

U = units2cgs;
P = physical_consts(U);
C = problem_specific_consts(U,P);
end

function P = physical_consts(U)

P.kB = 1.3806503e-16*U.erg/U.K; % Boltzmann const
P.NA = 6.022e12; % Avogadro's number
P.c = 2.99792458e10*U.cm/U.s; % speed of light
P.hP = 6.62606876e-27*U.erg*U.Hz; P.h = P.hP; % Plank's constant
P.hbar = P.hP/(2*pi); P.hb = P.hbar; % reduced Plank's constant
P.re = 2.816e-13*U.cm; % classical elecctron radius
P.e = 1.602176462e-19*U.C; % electron charge
P.mu = 1.66e-24*U.gr; % unified atomic mass unit
P.me = 9.10938188e-31*U.kg; % electron mass
P.mp = 1.66053873e-27*U.kg; % proton mass
P.gs = 2.002; % electron g-factor
P.uB = P.e*P.hbar/(2*P.me*P.c); % Bohr magneton
P.uN = P.e*P.hbar/(2*P.mp*P.c); % Nuclear magneton
P.ge = P.gs*P.uB/P.hbar; % electron gyromagnetic ratio (in [radians/(sec*gauss)])
P.gN = -3.826*P.uN/P.hbar; % neutron gyromagnetic ratio (in [radians/(sec*gauss)])
% P.eps0 = 8.854e-12*U.F/U.m; % MKS vacuum permitivity
% % P.eps0 = 1e7/(4*pi)*(U.Amp^2/U.N)/P.c^2; % MKS vacuum permitivity
% P.mu0 = 4*pi*1e-7*(U.N/U.Amp^2); % MKS vacuum susptibility
P.eps0 = 1/(4*pi); % CGS vacuum permitivity
P.mu0 = (4*pi)/P.c^2; % CGS vacuum susptibility
P.a0 = 4*pi*P.eps0*P.hbar^2/(P.e^2*P.me); % Bohr radius
P.re = (P.e)^2/(P.me*(P.c)^2);
P.eta0 = sqrt(P.mu0/P.eps0);
end

function C = problem_specific_consts(U,P)

C.gRb85 = P.ge/(2*5/2+1) ; % gyromagnetic ratio (in [radians/(sec*gauss)])
C.gRb87 = P.ge/(2*3/2+1) ; % gyromagnetic ratio (in [radians/(sec*gauss)])
C.g1H = 42.576*1E6*(2*pi)*U.rad/(U.s*U.T);
C.I1H = 1/2;

% C.g3He = -32.434*1e6*(2*pi)*U.rad/(U.s*U.T);
C.u3He = -2.1277*P.uN; % -1.0639*P.uN;
% C.g3He = C.u3He / (1/2);
C.g3He = -0.7618*C.g1H;
C.I3He = 1/2;

% C.u129Xe = 0.7729*P.uN;
C.u129Xe = -0.778*P.uN;
% C.g129Xe = C.u129Xe / (1/2)
C.g129Xe = -0.2781*C.g1H;
C.I129Xe = 1/2;

% C.g131Xe = +32.8E6/9.4E4*2*pi*(1.184/1.175);
% C.I131Xe = 3/2;

C.u131Xe = +0.692*P.uN;
% C.g131Xe = C.u131Xe / (3/2);
C.g131Xe = +0.0824*C.g1H;
C.I131Xe = 3/2;

C.u21Ne = -0.6617*P.uN;

C.dens_temp_A_K = 4.402; % from Seltzer thesis
C.dens_temp_B_K = 4453; % from Seltzer thesis
C.dens_temp_A_Rb = 4.312; % from Seltzer thesis
C.dens_temp_B_Rb = 4040; % from Seltzer thesis
C.dens_temp_A_Cs = 4.165; % from Seltzer thesis
C.dens_temp_B_Cs = 3830; % from Seltzer thesis

C.lambda_D1_K  = 770.1*U.nm; C.f_D1_K  = 1/3; % according to Kornak
C.lambda_D2_K  = 766.7*U.nm; C.f_D2_K  = 2/3; % according to Kornak
C.lambda_D1_Rb = 795.0*U.nm; C.f_D1_Rb = 0.342; % according to Steck
C.lambda_D2_Rb = 780.2*U.nm; C.f_D2_Rb = 0.696; % according to Steck

C.HF_A_K  = 2*pi*0.2309*U.GHz; C.HF_G_K  = 2*C.HF_A_K ;
C.HF_A_Rb = 2*pi*6.3417*U.GHz; C.HF_G_Rb = 2*C.HF_A_Rb;

end
