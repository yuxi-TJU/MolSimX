    function IH=HoppingR(betah,Vh)
    %hbar=1.055e-34;
    q=1.602e-19;
    %I0=q*q/hbar;
    kB=0.000086;
    Eg=betah(1);
    lambda=betah(2);
    T=298;              %×óµç¼«¾àÀë
    gr=betah(3);
    gl=betah(4);
    %rou=10;
    %gr=((2*pi)/hbar)*gamar^2*rou;
    %gl=((2*pi)/hbar)*gamal^2*rou;
    %c=beta(5);
    alphad1=betah(5);%%%×ó²à·ÖÑ¹
    alphad=1-alphad1;
    %Energy grid
    NE=16001;
    E=linspace(-8,8,NE);
    dE=E(2)-E(1);
    IV=length(Vh);
    % FC factor in Marcus transfer rate equation
    %FCab=(1/(4*pi*lambda*kB*T)^0.5).*exp(-(E-eg-lambda).^2./(4*lambda*kB*T));
    %FCab=FCab./(dE*sum(FCab));
    %FCba=(1/(4*pi*lambda*kB*T)^0.5).*exp(-(-E+eg-lambda).^2./(4*lambda*kB*T));
    %FCba=FCba./(dE*sum(FCba));
    %Bias
    for iV=1:IV
        FCrab=(1/(4*pi*lambda*kB*T)^0.5).*exp(-(lambda+((Eg-(alphad)*Vh(iV))-E)).^2./(4*lambda*kB*T));
        FCrba=(1/(4*pi*lambda*kB*T)^0.5).*exp(-(lambda-((Eg-(alphad)*Vh(iV))-E)).^2./(4*lambda*kB*T));
        FClab=(1/(4*pi*lambda*kB*T)^0.5).*exp(-(lambda+((Eg+(1-alphad)*Vh(iV))-E)).^2./(4*lambda*kB*T));
        FClba=(1/(4*pi*lambda*kB*T)^0.5).*exp(-(lambda-((Eg+(1-alphad)*Vh(iV))-E)).^2./(4*lambda*kB*T));
        flab=1./(1+exp((E)./(kB*T)));
        flba=1./(1+exp(-(E)./(kB*T)));
        frba=1./(1+exp(-(E)./(kB*T)));
        frab=1./(1+exp((E)./(kB*T)));
        %Rate integration
        Rabl=dE*sum(flab.*FClab)*gl;
        Rbal=dE*(sum(flba.*FClba))*gl;
        Rabr=dE*sum((frab.*FCrab))*gr;
        Rbar=dE*(sum(frba.*FCrba))*gr;
        IH(iV)=20000*(Rabl*Rbar-Rabr*Rbal)/(Rabl+Rbar+Rabr+Rbal)*q;% Current
    end
    end