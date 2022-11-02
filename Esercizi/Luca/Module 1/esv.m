a=2.4650e+4;
mu=398600;
e=0.7108;
t0=0;
t=2*pi*sqrt(a^3/mu);
f0=3*pi;
e=0.7108;
E0=2*atan(sqrt((1-e)/(1+e))*tan(f0/2));
if (E0<0)
    E0=2*pi-abs(E0);
end
M0= E0-e*sin(E0);
M= M0+sqrt(mu/a^3)*(t-t0);
if (M>0)
k=floor(M/(2*pi));
else 
k=abs(ceil(M/(2*pi)));
end
Mseg=wrapTo2Pi(M);
f=@(E) M-E+e*sin(E);
Eg= Mseg+(e*sin(Mseg)/(1-sin(Mseg+e)+sin(Mseg)));
options=optimset('TolX',1e-8);
E=fzero(f,Eg,options);
fseg=acos((cos(E)-e)/(1-e*cos(E)));
f=fseg+2*k*pi;