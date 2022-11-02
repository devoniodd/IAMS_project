function df= ode_orbit2BP_perturbed(~,f,mu,J2,R)
x=f(1);
y=f(2);
z=f(3);
vx=f(4);
vy=f(5);
vz=f(6);
r=norm([x y z]);
ax=-mu*x/r^3+((3/2)*J2*mu*R^2/r^4)*(x/r*((5*z^2/r^2)-1));
ay=-mu*y/r^3+((3/2)*J2*mu*R^2/r^4)*(y/r*((5*z^2/r^2)-1));
az=-mu*z/r^3+((3/2)*J2*mu*R^2/r^4)*(z/r*((5*z^2/r^2)-3));
df=[vx;vy;vz;ax;ay;az];
end