%% Orbit-propagator test
clear
close all
clc

load("../Data/utils.mat",'mu');

O = [42177.39,0.2,deg2rad(15),pi/4,deg2rad(80),0,pi; 
   2.9806e+04,0.6980,deg2rad(15),pi/4,deg2rad(80),pi,2*pi;
    9000,0.0001,deg2rad(15),pi/4,deg2rad(80),0,deg2rad(100);
    9000,0.0001,pi/4,pi/4,deg2rad(80),deg2rad(100),2*pi;
    9000,0.0001,pi/4,pi/4,deg2rad(80),0,deg2rad(30);
    19500,0.5385,pi/4,pi/4,deg2rad(120),0,pi
    ];

O2 = [50000.68854055203  0.80865268303717344  0.461525983980154  1.17569878936122  0.982554186272932  0  2*pi];

T = 2*pi/sqrt(mu) * O2(1)^(3/2)

t_step = 50;
th_step = 0.01;
warpfactor = 9000;
[lat,lon,CAR] = orbitpropagator(O2,t_step,warpfactor);
