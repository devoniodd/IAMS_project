function df= ode_orbit2BP(~,f,mu)
x=f(1);
y=f(2);
z=f(3);
vx=f(4);
vy=f(5);
vz=f(6);
r=norm([x y z]);
ax=-mu*x/r^3;
ay=-mu*y/r^3;
az=-mu*z/r^3;
df=[vx;vy;vz;ax;ay;az];
end