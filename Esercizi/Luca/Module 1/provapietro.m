a=7000;
mu=398600;
e=0.7108;
t0=0;
P=(2.000001)*pi*sqrt(a^3/mu);
f0=0;
e=[0,0.2,0.4,0.4,0.6,0.8,0.95];
t=linspace(t0,2*P,200);
for j=1:length(e)
for i=1:length(t)
y(i)=trueanomaly(a,e(j),mu,t0,f0,t(i));
end
f(:,j)= y;
end

figure(2)
[X,Y]=meshgrid(e,t);
Z=f;
surf(X,Y,Z)

