function t=Be10Newton(Te,x,D,n,attenuation,density,decay)
% Newton methods to find t using eroded thickness (D). 
% Te: effective age
% x: first guess of the age
% D: eroded thickness
% n: number of iterations
% attenuation: attenuation length of nucleon
% density: sample sediment density
% decay: decay constant
B=-density*D/attenuation;
f=@(x)exp(B)*exp(-decay*x)-B*Te/x+(decay*Te-1);
df=@(x)-decay*exp(B-decay*x)+B*Te/(x^2);
for i=1:n
    x=x-f(x)/df(x);
end
t=x;