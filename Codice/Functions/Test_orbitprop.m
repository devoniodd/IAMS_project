%% Orbit-propagator test
clear
close all
clc

%O = [42177.39,0.2,deg2rad(15),pi/4,deg2rad(80),0,pi; 
%   2.9806e+04,0.6980,deg2rad(15),pi/4,deg2rad(80),pi,2*pi;
%    9000,0.0001,deg2rad(15),pi/4,deg2rad(80),0,deg2rad(100);
%    9000,0.0001,pi/4,pi/4,deg2rad(80),deg2rad(100),2*pi;
%    9000,0.0001,pi/4,pi/4,deg2rad(80),0,deg2rad(30);
%    19500,0.5385,pi/4,pi/4,deg2rad(120),0,pi
%    ];
% O = [9723.68854055203  0.0865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  1.40668276722613  6.28318530717959;
%    -4061.20207150079  1.63674271526380  0.334061206030122  1.64741510007590  3.69717622068166  4.53830038309995  6.28318530717959;
%    16720  0.250200000000000  1.11900000000000  0.624500000000000  3.13500000000000  3.10000000000000  6.28318530717959];
O = 


t_step = 50;
warpfactor = 9000;
[lat,lon,CAR] = orbitpropagator(O,t_step,warpfactor);
