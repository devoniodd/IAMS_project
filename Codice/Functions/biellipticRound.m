function [dV,tT,tOrbit] = biellipticRound(ri,rf,rMax)

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% INITIAL AND TARGET ORBITS

Vi = sqrt(mu/ri);
Vf = sqrt(mu/rf);

%% TRANSFER ORBIT 1 (first ellipse)

hT1 = sqrt(2*mu)*sqrt(ri*rMax/(ri+rMax));
aT1 = (ri+rMax)/2;

VT1p = hT1/ri;
VT1a = hT1/rMax;

tT1 = pi * sqrt((aT1^3)/mu);

%% TRANSFER ORBIT 2 (second ellipse)

hT2 = sqrt(2*mu)*sqrt(rMax*rf/(rf+rMax));
aT2 = (rf+rMax)/2;

VT2p = hT2/rf;
VT2a = hT2/rMax;

tT2 = pi * sqrt((aT2^3)/mu);

%% OUTPUTS

dV = [VT1p-Vi; VT2a-VT1a; Vf-VT2p];

tT = [tT1; tT2];

tOrbit1.e = (rMax-ri)/(rMax+ri);
tOrbit1.a = aT1;
tOrbit1.p = aT1*(1-tOrbit1.e^2);
tOrbit1.ra = rMax;
tOrbit1.rp = ri;
tOrbit1.Vp = VT1p;
tOrbit1.Va = VT1a;

tOrbit2.e = (rMax-rf)/(rMax+rf);
tOrbit2.a = aT2;
tOrbit2.p = aT2*(1-tOrbit2.e^2);
tOrbit2.ra = rMax;
tOrbit2.rp = rf;
tOrbit2.Vp = VT2p;
tOrbit2.Va = VT2a;

tOrbit = [tOrbit1; tOrbit2];



