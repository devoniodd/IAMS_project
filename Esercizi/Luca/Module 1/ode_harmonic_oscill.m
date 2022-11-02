function dy= ode_harmonic_oscill(~,y,w0,gamma)
dy= [y(2); -gamma*y(2)-w0*y(1)]
end