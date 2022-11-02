function drawOrbit(X,Y,Z,V,timescale)
% Description to be added
% X,Y,Z position vectors
% V speed vector

%% Creating a time reference system
for i = 1:step_animation:length(X)-1
    ds = norm([X(i),Y(i),Z(i)] - [X(i+step_animation),Y(i+step_animation),Z(i+step_animation)]);
    dt(i) = (ds/V(i))/timeScaleFactor;
end