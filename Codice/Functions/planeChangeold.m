function [dV,targetOrbit] = planeChangeold(currentOrbit,targetOrbit,pureChange)

%% UTILS
if ismac
    load("../Data/utils.mat",'mu');
else
    load("..\Data\utils.mat",'mu');
end

%% INPUT ORGANIZATION
ai = currentOrbit(1);
ei = currentOrbit(2);
ii = currentOrbit(3);
Oi = currentOrbit(4);
oi = currentOrbit(5);

af = targetOrbit(1);
ef = targetOrbit(2);
iF = targetOrbit(3);
Of = targetOrbit(4);
if ~pureChange
    of = targetOrbit(5);
end

%% PROBLEM VARIABLES

dO = Of - Oi;
di = iF - ii;

%% CASES 

if dO > 0 && di > 0
    id = 1;
elseif dO > 0 && di < 0
    id = 2;
elseif dO < 0 && di > 0
    id = 3;
elseif dO < 0 && di < 0
    id = 4;
else 
    error("Error, dO and/or di are not valid, please check inputs!")
end

%% SPHERIC TRIANGLES
alpha = acos(cos(ii)*cos(iF) + sin(ii)*sin(iF)*cos(dO));
uiCos= (-1)^(id+1) * (cos(alpha)*cos(ii) - cos(iF))/(sin(alpha)*sin(ii));
uiSin = (-1)^(id+1) * sin(iF)*sin(dO)/sin(alpha);
ui = atan(uiSin/uiCos);
ufCos = (-1)^(id+1) * (cos(ii) - cos(alpha)*cos(iF))/(sin(alpha)*sin(iF));
ufSin = (-1)^(id+1) * sin(ii)*sin(dO)/sin(alpha);
uf = atan(ufSin/ufCos);

if ui >= oi
    thetai = ui - oi;
else 
    thetai = 2*pi + ui - oi;
end

if pureChange
    thetaf = thetai;
  if uf >= thetaf
        of = uf - thetaf;
    else 
        of = 2*pi + uf - thetaf;
    end
else
    if uf >= of
        thetaf = uf - of;
    else 
        thetaf = 2*pi + uf - of;
    end
end

Vit = sqrt(mu/(ai*(1-ei^2))) * (1 + ei*cos(thetai));

if pureChange
    dV = 2 * Vit * sin(alpha/2);
    return;
end

Vir = sqrt(mu/(ai*(1-ei^2))) * ei * sin(thetai);

Vft = sqrt(mu/(af*(1-ef^2))) * (1 + ef*cos(thetaf));
Vfr = sqrt(mu/(af*(1-ef^2))) * ef * sin(thetaf);

dV = sqrt((Vfr-Vir)^2 + Vit^2 + Vft^2 - 2*Vit*Vft*cos(alpha));