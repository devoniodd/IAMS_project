function [dV,tT] = biellipticRound(ri,rf,rMax)

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



