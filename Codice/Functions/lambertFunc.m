function [v1,v2,z] = lambertFunc(r1,r2,dt,prograde)

%% UTILS IMPORT 
addpath(genpath("../Data"));
addpath(genpath("Stumpff"));
load("utils.mat",'mu','kDir');

r1Vec = r1;
r2Vec = r2;
r1 = norm(r1);
r2 = norm(r2);

%% dTH

discern = dot(cross(r1Vec,r2Vec),kDir);
if prograde == 1
    if discern > 0
        dTh = acos(dot(r1Vec,r2Vec) / (r1 * r2));
    else
        dTh = 2*pi - acos(dot(r1Vec,r2Vec) / (r1 * r2));
    end
else
    if discern < 0
        dTh = acos(dot(r1Vec,r2Vec) / (r1 * r2));
    else
        dTh = 2*pi - acos(dot(r1Vec,r2Vec) / (r1 * r2));
    end
end

%% A

A = sin(dTh) * sqrt((r1*r2)/(1-cos(dTh)));

y = @(z) r1 + r2 + A*(z*stumpS(z) - 1)/sqrt(stumpC(z));

F = @(z) (y(z)/stumpC(z))^1.5 * stumpS(z) + A*sqrt(y(z)) - sqrt(mu) * dt;

guess = 20;
while F(guess) > 0
    guess = guess - 5;
end

while F(guess) < 0
    guess = guess + 0.1;
end
optnew = optimset('TolX',1e-15);
z = fzero(F,guess,optnew);

f = 1-y(z)/r1;
g = A*sqrt(y(z)/mu);
dg = 1 - y(z)/r2;
v1 = 1/g * (r2Vec - f*r1Vec);
v2 = 1/g * (dg * r2Vec - r1Vec);

        

end

