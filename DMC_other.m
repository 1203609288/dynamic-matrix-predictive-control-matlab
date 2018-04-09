clear all;
clc;
close all;
d=35;
A=[1.6375 -0.6703];
num=[0.0350 0.0307];
B=[zeros(1,d) num];
den=[1 -A zeros(1,d-1)];
P=tf(num,den,0.2)
Ap=A;Bp=B;
na=size(Ap,2);
nb=size(Bp,2);
y_ant=zeros(1,na);
u_ant=zeros(1,nb);
N=80;               %模型长度     
g=step(P,N);
dulib=zeros(1,N);
p=60;               %预测时域
m=3;                %控制时域
G=zeros(p,m);
G(:,1)=g(1:p);
for i=2:m
    ga=[zeros(1,i-1) g(1:p-i+1)'];
    G(:,i)=ga';
end
%--------define simulation---------------
puntos=1000;         %仿真结束时间
lambda=1;           %控制量加权
Mn=(G'*G+lambda*eye(m));
Mn1=inv(Mn)*G';
qn=Mn1(1,:);
u=zeros(1,puntos);
yp=zeros(1,puntos);
ypr=zeros(1,puntos);
w=[zeros(1,5) 0.5*ones(1,puntos-5)];
af=0.1;            %输入滤波系数,越小性能越好
r=filter((1-af),[1 -af],w);
%--------inicialization-------------------
inc_u=0;
for k=1:puntos
    yp(k)=Ap*y_ant'+Bp*u_ant';
    %----add a disturbance at the output
    if k>100
        ypr(k)=yp(k)-0.1;
    else
        ypr(k)=yp(k);
    end
    %----actualization of the output
    if na==1
        y_ant=yp(k);
    else
        aux_y=y_ant(1:na-1);
        y_ant=[yp(k) aux_y];
    end
    %----computes prediction for DMC
    f=zeros(1,p);
    for kk=1:p
        for i=1:N-p
            vect_g(i)=g(kk+i)-g(i);
        end
        for i=N-p+1:N
            vect_g(i)=g(N)-g(i);
        end
        f(kk)=ypr(k)+vect_g*dulib';
    end
    ref=r(k)*ones(1,p)';
    inc_u=qn*(ref-f');
    if k==1
        u(k)=inc_u;
    else
        u(k)=u(k-1)+inc_u;
    end
    %actualization for control vector
    aux_u=u_ant(1:nb-1);
    u_ant=[u(k) aux_u];
    %actualization dulib
    aux_2=dulib(1:N-1);
    dulib=[inc_u aux_2];
end

%------plot-------
nm=puntos;
h=subplot(2,1,1);
for i=1:nm
    t(i)=0.2*i;
end
plot(t(1:nm),w(1:nm),'-.',t(1:nm),ypr(1:nm),'-','LineWidth',1.5);
title('Output of the process');
legend('r','DMC','Location','SouthEast');
xlabel('time(s)','FontSize',12);
ylabel('r,y','FontSize',12);
axis([0 50 -0.1 0.7]);
set(h,'FontSize',12);
grid on;
h=subplot(2,1,2);
plot(t(1:nm),u(1:nm),'-','LineWidth',1.5);
title('Control Action of the process');
xlabel('time(s)','FontSize',12);
ylabel('u','FontSize',12);
axis([0 50 -0.5 0.5]);
set(h,'FontSize',12);
grid on;



        