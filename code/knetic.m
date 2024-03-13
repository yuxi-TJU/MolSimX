function [result] = knetic(VV,scanr,kn,kp,oxi1,red1,oxi2,red2,I1,I2)
q=1.602e-19;
kB=0.000086;

T=300; %�¶�
poi = 100*VV+1; %һ����ѹȡ����
Vh1=linspace(0,VV,poi);%��һ��ɨ��ĵ�ѹ���飨0V-2V��
Vh2=linspace(VV,0,poi);%�ڶ���ɨ��ĵ�ѹ���飨2V-0V��
Vh3=linspace(0,-VV,poi);%������ɨ��ĵ�ѹ���飨0V--2V��
Vh4=linspace(-VV,0,poi);%������ɨ��ĵ�ѹ���飨-2V-0V��
Vhall=[Vh1,Vh2,Vh3,Vh4];%���в���d ���飨0V-2V--2V-0V��
dt=(VV/(poi-1))/scanr; %ʱ�����
t=0;%��ʼʱ��
AQ=1;%��������
H2AQ=0;%���ӻ���������
NE=501;
E=linspace(-5,5,NE);
dE=E(2)-E(1);
IV=2*poi;

%��һ����0V-2V���еı仯
for i = 1:poi
    m1 = oxi1(poi+i);%����̬
    m2 = red1(poi+i);%��ԭ̬
    n1 = oxi2(poi+i);%���ӻ�����̬
    n2 = red2(poi+i);%���ӻ���ԭ̬
    
    t = t + dt;%ʱ��ı仯��ʱ��仯dA�仯��
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
    
    m1 = oxi1(2*poi+1-i);%����̬
    m2 = red1(2*poi+1-i);%��ԭ̬
    n1 = oxi2(2*poi+1-i);%���ӻ�����̬
    n2 = red2(2*poi+1-i);%���ӻ���ԭ̬
    t = t + dt;%ʱ��ı仯��ʱ��仯dA�仯��
    kn1=kn*n1;
    kp1=kp*m2;
    AQ=(AQ-(kn1/(kp1+kn1)))*exp(-(kp1+kn1)*dt)+(kn1/(kp1+kn1));
    H2AQ=1-AQ;
    amoAQ2(i) = AQ;
    amoh2AQ2(i) = H2AQ;
    CurrentsAQ2(i) = I1(2*poi+1-i)*AQ;
    Currentsh2AQ2(i) = I2(2*poi+1-i)*H2AQ;
end
%��������0V--2V���еı仯
for i = 1:poi
    
    m1 = oxi1(poi+1-i);%����̬
    m2 = red1(poi+1-i);%��ԭ̬
    n1 = oxi2(poi+1-i);%���ӻ�����̬
    n2 = red2(poi+1-i);%���ӻ���ԭ̬
    t = t + dt;%ʱ��ı仯��ʱ��仯dA�仯��
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
    m1 = oxi1(i);%����̬
    m2 = red1(i);%��ԭ̬
    n1 = oxi2(i);%���ӻ�����̬
    n2 = red2(i);%���ӻ���ԭ̬
    t = t + dt;%ʱ��ı仯��ʱ��仯dA�仯��  
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
