function S = stumpS(z)

if z > 0
    S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
elseif z < 0 
    S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
else
    S = 1/6;
end

end