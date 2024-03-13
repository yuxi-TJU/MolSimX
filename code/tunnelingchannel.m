function[result1] = tunnelingchannel(Eg,u, g1,g2,V,alpha,T)
e=1.60217*(10^(-19));
h=6.62607004*(10^(-34));
hbar=h/(2*pi);

kB=0.000086;
kT=kB*T;
g=g1+g2;


NE=16001;
E=linspace(-8,8,NE);
dE=E(2)-E(1);

D=(g/(2*pi))./(((E-(Eg+u)).^2)+((g/2)^2));% Lorentzian Density of states per eV
D=D./(dE*sum(D));%Normalizing to one
UL=(alpha*V);
UR=((1-alpha)*V);
f1=1./(1+exp((E+UL)./kT));
f2=1./(1+exp((E-UR)./kT));
Current=-dE*e*e/h*(sum(D.*(f1-f2)))*2*pi*(g1*g2/g);
N=dE*sum(D.*((g1*f1+g2*f2)/g));
ox=1-N;
result1=[V,Current,N,ox];
end
