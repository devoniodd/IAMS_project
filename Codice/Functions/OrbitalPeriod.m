function [T] = orbitalPeriod(a,M)
% The function allows to calculate the orbital period T from the radius.
% INPUT:
% a - ellipse major apsis
% M - mass of the planet (If left unspecified erth mass is used)
%
% OUTPUT:
% T - orbital period

% To call this function type[ P ] = orbital period ( R )

%% UTILS IMPORT

if ismac
    load("../Data/utils.mat",'G','mEarth');
else
    load("..\Data\utils.mat",'G','mEarth');
end

%% PERIOD

if nargin == 1
    M = mEarth;
end

T = sqrt((4*(pi^2)*a^3)/(G*M));

end