clear all;
close all;
clc;

r1 = [-8048.2861, -4171.3048, 2895.9296];
v1 = [1.7540, -5.9910, -1.9520];

orbit1 = zeros(1,7);
[orbit1(1), e, orbit1(3), orbit1(4), orbit1(5), orbit1(6)] = carToOrbital(r1,v1);
orbit1(2) = norm(e);

orbit2 = [16720.0000 0.2502 1.1190 0.6245 3.1350 3.1000];
[r2, v2] = orbitalToCar(16720.0000,0.2502,1.1190,0.6245,3.1350,3.1000);
r2 = r2';

v1t = [30, 30, 30];
t = 0;
orbitT = zeros(1,7);
orbitT(1) = 1000;
orbitT(2) = 0.5;
vLim = 5;
while orbitT(1) * (1-orbitT(2)) < 6578
    t = t+200;
    while (abs(norm(v1) - norm(v1t))) > vLim
        [v1t, v2t] = lambert(r1,r2,t,1);
        t = t + 10;
    end
    orbitT = zeros(1,7);
    [orbitT(1), e, orbitT(3), orbitT(4), orbitT(5), orbitT(6)] = carToOrbital(r1,v1t);
    orbitT(2) = norm(e);
    v1t = [30, 30, 30];
end

orbit1(7) = 2*pi;
orbit2(7) = 2*pi;
orbit2(6) = 0;
orbitT(7) = pi;

orbits = [orbit1; orbitT; orbit2];
