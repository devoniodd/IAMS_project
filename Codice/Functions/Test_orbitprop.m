%% Orbit-propagator test
clear
close all
clc

O = [42177.39,0.2,deg2rad(15),pi/4,deg2rad(80),0,2*pi];
t_step = 50;
warpfactor = 9000;
[lon,TH] = orbitpropagator(O,t_step,warpfactor);
