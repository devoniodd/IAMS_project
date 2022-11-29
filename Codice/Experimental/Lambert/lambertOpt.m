clear all;
close all;
clc;

addpath(genpath("../../Data/"))
load("utils.mat",'mu');

%% ORBIT 1
r1 = [-8048.2861, -4171.3048, 2895.9296];
v1 = [1.7540, -5.9910, -1.9520];

orbit1 = zeros(1,7);
[orbit1(1), orbit1(2), orbit1(3), orbit1(4), orbit1(5), orbit1(6)] = carToOrbital(r1,v1);
e1 = orbit1(2);
a1 = orbit1(1);
p1 = a1 * (1-e1^2);
i1 = orbit1(3);
O1 = orbit1(4);
orb1Part = [orbit1(1),orbit1(2),orbit1(3),orbit1(4),orbit1(5),orbit1(6)];

%% ORBIT 2
orbit2 = [16720.0000 0.2502 1.1190 0.6245 3.1350 3.1000];
[r2, v2] = orbitalToCar(16720.0000,0.2502,1.1190,0.6245,3.1350,3.1000);
r2 = r2';

orb2Part = [orbit2(1),orbit2(2),orbit2(3),orbit2(4),orbit2(5),orbit2(6)];

%% FORs

t = 1000;
a = 4000;
e = 0.5;
vLim = 11;
orbitT = zeros(1,7);

it1 = 0;
it2 = 0;
n = round(2*pi/0.2);
dVs = struct();
ts = struct();
parfor thetaOpt1=0:1:round(2*pi * (1/0.2))+1
    thetaOpt1
    dVlambert = vLim;
    t = 1000;
    a = 4000;
    e = 0.5;
    [r1,v1] = orbitalToCar(9723.68854055203,0.086526830371734,0.461525983980154,1.17569878936122,0.982554186272932,thetaOpt1/5);
    v1t = sqrt(mu/p1) * (1+e1*cos(thetaOpt1/5));
    v1r = sqrt(mu/p1) * (e1*sin(thetaOpt1/5));
    it2 = 0;
    for thetaOpt2=0:0.2:2*pi
        it2 = it2 + 1
        [r2,v2] = orbitalToCar(16720,0.250,1.119,0.6245,3.13500,thetaOpt2);
        itIns = 0;
        while a * (1-e) < 6578
            itIns = itIns + 1; 
            if itIns < 100
                t = t+10;
            elseif itIns < 500
                t = t+100;
            else
                t = t+200;
            end
            vt1= lambert(r1,r2,t,1);
            [a, e, i, O] = carToOrbital(r1,vt1);
        end
        itIns = 0;
%         while dVlambert >= vLim
%             itIns = itIns + 1; 
%             if itIns < 100
%                 t = t+10;
%             elseif itIns < 500
%                 t = t+100;
%             else
%                 t = t+200;
%             end
%             [vt1, vt2] = lambert(r1,r2,t,1);
%             [at, et, it, Ot, ~, thetat] = carToOrbital(r1,vt1);
%             pt = at * (1-et^2);
%             vtt = sqrt(mu/pt) * (1+et*cos(thetat));
%             vtr = sqrt(mu/pt) * (et*sin(thetat));
%             alpha = acos(cos(i1)*cos(it) + sin(i1)*sin(it)*cos(O1-Ot));
%             dVlambert = sqrt((vtr-v1r)^2+v1t^2+vtt^2-2*vtt*v1t*cos(alpha));
%         end
        [at, et, it, Ot, ~, thetat] = carToOrbital(r1,vt1);
            pt = at * (1-et^2);
            vtt = sqrt(mu/pt) * (1+et*cos(thetat));
            vtr = sqrt(mu/pt) * (et*sin(thetat));
            alpha = acos(cos(i1)*cos(it) + sin(i1)*sin(it)*cos(O1-Ot));
            dVlambert = sqrt((vtr-v1r)^2+v1t^2+vtt^2-2*vtt*v1t*cos(alpha));
        dVs(thetaOpt1+1).matrix(thetaOpt1,it2) = dVlambert;
        ts(thetaOpt1+1).matrix(thetaOpt1,it2) = t;
        a = 4000;
        e = 0.5;
        orbitT = zeros(1,7);
        t = 1000;
        dVlambert = vLim;
    end

end

