function IT=TunnelingR(beta,V)
    e=1.60217*(10^(-19));
    h=6.62607004*(10^(-34));
    hbar=h/(2*pi);
    I0=e*e/hbar;
    %T=298;
    %g1=0.01;
    IV=length(V);
    Eg=beta(1);
        g1=beta(2);
        % c=beta(3);
        g2=beta(3);
        %c=beta(4);
        g=g1+g2;
        %c=g1*g2/(g1+g2);
        alpha=beta(4);
        %½âÎö½â
        NE=1601;
        E=linspace(-8,8,NE);
        dE=E(2)-E(1);
        D=(g/(2*pi))./((E.^2)+((g/2)^2));% Lorentzian Density of states per eV
        D=D./(dE*sum(D));%Normalizing to one
        %Bias
        %kT=(k*T)/e;
        kT=0.025;
        for iV=1:IV
            Vd=V(iV);
            UL=(alpha*Vd);
            UR=((1-alpha)*Vd);
            f1=1./(1+exp((E+Eg+UL)./kT));
            f2=1./(1+exp((E+Eg-UR)./kT));
            IT(iV)=2000*dE*I0*(sum(D.*(f1-f2)))*(g1*g2/g);
        end
    end
