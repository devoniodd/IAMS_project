a=2.4650e+4;
mu=398600;
e=0.7108;
t0=0;
t=2*pi*sqrt(a^3/mu);
f0=3*pi;
e=0.7108;
E0=acos((cos(f0)+e)/(1+e*cos(f0)));
if(f0>=2*pi)
k=floor(f0/(2*pi));
E0=E0+2*k*pi;
end
M0= E0-e*sin(E0);
M= M0+sqrt(mu/a^3)*(t-t0);
Mseg=wrapTo2Pi(M);
if (M>0)
j=floor(M/(2*pi));
else 
j=abs(ceil(M/(2*pi)));
end
M=Mseg+2*j*pi;
f=@(E) Mseg-E+e*sin(E);
Eg= Mseg+(e*sin(Mseg)/(1-sin(Mseg+e)+sin(Mseg)));
options=optimset('TolX',1e-8);
E=fzero(f,Eg,options);
fseg=2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
f=fseg+2*j*pi;