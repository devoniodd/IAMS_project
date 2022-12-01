clear all;
close all;
clc;

%% INPUTS

addpath(genpath("../../Data/"))
load("utils.mat",'mu');
load("PForbs.mat");

%% INITIAL AND FINAL ORBIT DEFINITION

[r1, V1] = orbitalToCar(orb1(1),orb1(2),orb1(3),orb1(4),orb1(5),orb1(6));

[r2, V2] = orbitalToCar(orb2(1),orb2(2),orb2(3),orb2(4),orb2(5),orb2(6));

%% LAMBERT PARAMETERS


%Parallelization parameters
delete(gcp('nocreate'));
parpool(4);

theta1 = 0:1:2*pi*10;
theta2 = 0:0.1:2*pi;

VVS = zeros(length(theta1));
TTS = zeros(length(theta1));

parfor th1 = theta1

    [r1, V1] = orbitalToCar(orb1(1),orb1(2),orb1(3),orb1(4),orb1(5),th1/10);

    % Starting time of flight
    t = 3000;

    % Exit parameters
    hMin = 200;

    % Utils
    r0 = 6378;
    orbitT = zeros(1,7);
    
    dVtotList = zeros(1,length(theta1));
    ts = zeros(1,length(theta1));


    for th2 = theta2
        
        it = 1;

        [r2, V2] = orbitalToCar(orb2(1),orb2(2),orb2(3),orb2(4),orb2(5),th2);
        while it == 1 || a*(1-e) < r0+hMin

            % Transfer orbit velocity
            [Vt1,Vt2] = lambertFunc(r1,r2,t,1);

            % Find needed orbital parameters
            [a, e] = carToOrbital(r1,Vt1);

            % Update
            it = it+1;
            t = t+10;

        end

        dVlambert = norm(Vt1-V1);
        dVfinal = norm(Vt2-V2);
        dVtot = dVlambert + dVfinal;

        it = 1;

        dVtotList(it) = norm(dVtot);
        ts(it) = t;

        

        while it == 1 || dVtotList(it) < dVtotList(it-1)

            % Update
            t = t+10;
            it = it+1;

            % Transfer orbit velocity
            [Vt1,Vt2] = lambertFunc(r1,r2,t,1);

            % Velocities
            dVlambert = norm(Vt1-V1);
            dVfinal = norm(Vt2-V2);
            dVtot = dVlambert + dVfinal;

            dVtotList(it) = norm(dVtot);
            ts(it) = t;

        end

    end

    VVS(th1,:) = dVtotList;
    TTS(th1,:) = ts;

end

delete(gcp('nocreate'));


% [orbitT(1), orbitT(2), orbitT(3), orbitT(4), orbitT(5), orbitT(6)] = carToOrbital(r1,Vt1);
% [~, ~, ~, ~, ~, orbitT(7)] = carToOrbital(r2,Vt2);
% orbit1(7) = 2*pi;
% orbit2(7) = 2*pi;
% orbit2(6) = 3.1;
% orbits = [orbit1; orbitT; orbit2];