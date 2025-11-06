function U = units2cgs

U.rad = 1; U.mrad = 1E-3; U.urad = 1E-6; U.deg = 1/360; U.mdeg = 1e-3*U.deg;

U.m = 1e2; U.cm = 1; U.mm = 1e-3*U.m; U.um = 1e-6*U.m; U.u = U.um;
U.nm = 1e-9*U.m; U.A = 1e-10*U.m; U.pm = 1e-12*U.m; U.km = 1e3*U.m;
U.inch = 25.4*U.mm;

U.cm2 = U.cm^2; U.cm3 = U.cm^3;

U.K = 1;

U.kg = 1e3; U.gr = 1; U.g = U.gr;

U.s = 1; U.ms = 1e-3; U.us = 1e-6; U.ns = 1e-9; U.ps = 1e-12*U.s; U.fs = 1e-15*U.s; U.as = 1e-18*U.s;
U.min = 60*U.s; U.h = 60*U.min;
U.Hz = 1; U.kHz = 1e+3; U.MHz = 1e+6; U.GHz = 1e+9; U.THz = 1e+12;
          U.mHz = 1e-3; U.uHz = 1e-6; U.nHz = 1e-9;
U.rad_s   = U.rad/U.s; U.mrad_s   = U.mrad/U.s; U.urad_s   = U.urad/U.s;
U.rad_sec = U.rad/U.s; U.mrad_sec = U.mrad/U.s; U.urad_sec = U.urad/U.s;

U.m_s = U.m/U.s; U.km_s = U.km/U.s; U.mm_us = U.mm/U.us; U.m_ms = U.m/U.ms;

U.cm2_s = (U.cm)^2/U.s;

U.G = 1; U.T = 1e4*U.G; U.mG = 1e-3*U.G; U.uG = 1e-6*U.G; U.nG = 1e-9*U.G;
U.pG = 1e-12*U.G; U.kG = 1e3*U.G;
U.mT = 1e-3*U.T; U.uT = 1e-6*U.T; U.nT = 1e-9*U.T;
U.pT = 1e-12*U.T; U.fT = 1e-15*U.T;

U.dyn = 1; U.dyne = U.dyn;
U.N = 1e5*U.dyn;
U.dyn_cm2 = U.dyn/(U.cm^2); U.Pa = U.N/(U.m^2);
U.bar = 1e5*U.Pa; U.mbar = 1e-3*U.bar; U.ubar = 1e-6*U.bar; U.mubar = 1e-6*U.bar;
U.atm = 1.01325*U.bar;
U.torr = U.bar/760; U.tor = U.torr;
U.mtorr = 1e-3*U.torr; U.mtor = 1e-3*U.tor;
U.erg = 1; U.J = 1e7*U.erg; U.eV = 1.6022*U.J;
U.W = U.J/U.s; U.mW = 1e-3*U.W;
U.statC = 1*sqrt(U.erg*U.cm); U.stC = U.statC; U.C = 2.9979e9*U.statC;
U.statAmp = 1*sqrt(U.dyn); U.stAmp = U.statAmp; U.Amp = U.C/U.statC*U.statAmp;
U.statV = 1*sqrt(U.erg/U.cm); U.stV = U.statV; U.V = (1e8/U.C)*U.statV;
U.Ohm = 1e9/(U.C^2)*(U.s/U.cm); U.ohm = U.Ohm; U.F = 1e-9*(U.C^2)*(U.cm);

U.amg = 2.6868e19/U.cm3; % Loschmidt constant (1bar/kB*273.15K)

U.sqrtHz = 1; U.Hz_sqrtHz = U.Hz/U.sqrtHz;
U.deg_sqrth = U.deg/sqrt(U.h); U.mdeg_sqrth = U.mdeg/sqrt(U.h);
U.deg_h = U.deg/U.h; U.mdeg_h = U.mdeg/U.h;

clear U.C;
end