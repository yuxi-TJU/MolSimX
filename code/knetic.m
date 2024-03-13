function [result] = knetic(VV,scanr,kn,kp,oxi1,red1,oxi2,red2,I1,I2)
q=1.602e-19;
kB=0.000086;

T=300; %温度
poi = 100*VV+1; %一步电压取点数
Vh1=linspace(0,VV,poi);%第一步扫描的电压数组（0V-2V）
Vh2=linspace(VV,0,poi);%第二步扫描的电压数组（2V-0V）
Vh3=linspace(0,-VV,poi);%第三步扫描的电压数组（0V--2V）
Vh4=linspace(-VV,0,poi);%第三步扫描的电压数组（-2V-0V）
Vhall=[Vh1,Vh2,Vh3,Vh4];%所有步数d 数组（0V-2V--2V-0V）
dt=(VV/(poi-1))/scanr; %时间变量
t=0;%初始时间
AQ=1;%蒽醌总量
H2AQ=0;%质子化蒽醌总量
NE=501;
E=linspace(-5,5,NE);
dE=E(2)-E(1);
IV=2*poi;

%第一步（0V-2V）中的变化
for i = 1:poi
    m1 = oxi1(poi+i);%氧化态
    m2 = red1(poi+i);%还原态
    n1 = oxi2(poi+i);%质子化氧化态
    n2 = red2(poi+i);%质子化还原态
    
    t = t + dt;%时间的变化（时间变化dA变化）
    kn1=kn*n1;
    kp1=kp*m2;
    AQ=(AQ-(kn1/(kp1+kn1)))*exp(-(kp1+kn1)*dt)+(kn1/(kp1+kn1));
    H2AQ=1-AQ;
    amoAQ1(i) = AQ;
    amoh2AQ1(i) = H2AQ;
    CurrentsAQ1(i) = I1(poi+i)*AQ;
    Currentsh2AQ1(i) = I2(poi+i)*H2AQ;
    
end
for i = 1:poi
    
    m1 = oxi1(2*poi+1-i);%氧化态
    m2 = red1(2*poi+1-i);%还原态
    n1 = oxi2(2*poi+1-i);%质子化氧化态
    n2 = red2(2*poi+1-i);%质子化还原态
    t = t + dt;%时间的变化（时间变化dA变化）
    kn1=kn*n1;
    kp1=kp*m2;
    AQ=(AQ-(kn1/(kp1+kn1)))*exp(-(kp1+kn1)*dt)+(kn1/(kp1+kn1));
    H2AQ=1-AQ;
    amoAQ2(i) = AQ;
    amoh2AQ2(i) = H2AQ;
    CurrentsAQ2(i) = I1(2*poi+1-i)*AQ;
    Currentsh2AQ2(i) = I2(2*poi+1-i)*H2AQ;
end
%第三步（0V--2V）中的变化
for i = 1:poi
    
    m1 = oxi1(poi+1-i);%氧化态
    m2 = red1(poi+1-i);%还原态
    n1 = oxi2(poi+1-i);%质子化氧化态
    n2 = red2(poi+1-i);%质子化还原态
    t = t + dt;%时间的变化（时间变化dA变化）
    kn1=kn*n1;
    kp1=kp*m2;
    AQ=(AQ-(kn1/(kp1+kn1)))*exp(-(kp1+kn1)*dt)+(kn1/(kp1+kn1));
    H2AQ=1-AQ;
    amoAQ3(:,i) = AQ;
    amoh2AQ3(:,i) = H2AQ;  
    CurrentsAQ3(i) = I1(poi+1-i)*AQ;
    Currentsh2AQ3(i) = I2(poi+1-i)*H2AQ;
end
for i = 1:poi
    m1 = oxi1(i);%氧化态
    m2 = red1(i);%还原态
    n1 = oxi2(i);%质子化氧化态
    n2 = red2(i);%质子化还原态
    t = t + dt;%时间的变化（时间变化dA变化）  
    kn1=kn*n1;
    kp1=kp*m2;
    AQ=(AQ-(kn1/(kp1+kn1)))*exp(-(kp1+kn1)*dt)+(kn1/(kp1+kn1));
    H2AQ=1-AQ;
    amoAQ4(i) = AQ;
    amoh2AQ4(i) = H2AQ;
    CurrentsAQ4(i) = I1(i)*AQ;
    Currentsh2AQ4(i) = I2(i)*H2AQ;
end
amoAQ =[amoAQ1,amoAQ2,amoAQ3,amoAQ4]';
amoh2AQ =[amoh2AQ1,amoh2AQ2,amoh2AQ3,amoh2AQ4]';
CurrentsAQ=[CurrentsAQ1,CurrentsAQ2,CurrentsAQ3,CurrentsAQ4]';
Currentsh2AQ=[Currentsh2AQ1,Currentsh2AQ2,Currentsh2AQ3,Currentsh2AQ4]';
Currents=CurrentsAQ+Currentsh2AQ;
result=[Vhall',amoAQ,amoh2AQ,CurrentsAQ,Currentsh2AQ,Currents];
end
