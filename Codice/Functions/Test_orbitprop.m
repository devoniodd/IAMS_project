%% Orbit-propagator test
clear
close all
clc

O = [42177.39,0.2,deg2rad(15),pi/4,deg2rad(80),0,pi; 
    2.9806e+04,0.6980,deg2rad(15),pi/4,deg2rad(80),pi,2*pi;
    9000,0.0001,deg2rad(15),pi/4,deg2rad(80),0,deg2rad(100);
    9000,0.0001,pi/4,pi/4,deg2rad(80),deg2rad(100),2*pi;
    9000,0.0001,pi/4,pi/4,deg2rad(80),0,deg2rad(30);
    19500,0.5385,pi/4,pi/4,deg2rad(120),0,pi
    ];
t_step = 50;
warpfactor = 9000;
[lat,lon,CAR] = orbitpropagator(O,t_step,warpfactor);
