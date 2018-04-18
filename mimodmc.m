clear;clc;
num={[5],[3];[6],[9]};
den={[3 1 3],[1 2 5];[2 1 7],[2 3 6]};
sys=tf(num,den);%模型传递函数
g11=poly2tfd(num{1,1},den{1,1},0,0);
g12=poly2tfd(num{1,2},den{1,2},0,0);
g21=poly2tfd(num{2,1},den{2,1},0,0);
g22=poly2tfd(num{2,2},den{2,2},0,0);
delta=0.5;%采样时间
P=12;M=6;m=2;p=2;N=40; %   M,P,m,p分别为控制时域长度，预测时域长度，输入量个数，输出量个数，N为建模时域
ny=2;
tfinal=500;
mymodel=tfd2step(tfinal,delta,ny,g11,g12,g21,g22);%计算阶跃响应模型

%作图，画出模型的阶跃响应曲线
figure(1)
subplot(2,2,1);
step(num{1,1},den{1,1});
title('u1-y1阶跃响应');
xlabel('time');
subplot(2,2,2);
step(num{1,2},den{1,2});
title('u1-y2阶跃响应');
xlabel('time');
subplot(2,2,3);
step(num{2,1},den{2,1});
title('u2-y1阶跃响应');
xlabel('time');
subplot(2,2,4);
step(num{2,2},den{2,2});
title('u2-y2阶跃响应');
xlabel('time');

ywt=[3,1];%Q矩阵
uwt=[400,300];%R矩阵
alpha=[1,1];%H矩阵
r=[1;2];%设定值
tend=500;%结束时间

%计算反馈校正H矩阵
H=[];
for i=1:p
    h=alpha(1,i)*ones(N,1);
    H=blkdiag(H,h);
end

%计算移位矩阵S
for i=1:p
    for j=1:N-1
        S((i-1)*N+j,(i-1)*N+j+1)=1;%次对角线元素为1
    end
    S((i*N),(i*N))=1;%右下角元素为1
end

%输出设定值R矩阵
R=[];
for i=1:p
    r=r(i,1)*ones(P,1);
    R=[R;r];
end

y_Real=zeros(p,tend);%实际输出
e=zeros(p,tend);%误差
y=zeros(p,tend);%期望
U=zeros(m,tend);%输入

for i=1:p
    for j=1:N
        y_N((i-1)*N+j,1)=0;
        y_N0((i-1)*N+j,1)=0;
    end
    for j=1:P
        y_P0((i-1)*P+j,1)=0;
    end
end
deltY=[];
deltU=[];
y_Ncor=[];


[kmpc,A,L,A_N,a,Q]=DMC_martixD(mymodel,uwt,ywt,M,P,m,p,N);%计算D矩阵

for i=1:1:tend
    e(:,i)=y_Real(:,i)-y(:,i);
    y_Ncor(:,i)=y_N(:,i)+H*e(:,i);
    y_N0(:,i)=S*y_Ncor(:,i);
    for j=1:p
        y_P0((j-1)*P+1:j*P,i)=y_N0((j-1)*N+1:(j-1)*N+P,i);
    end
    deltY(:,i)=R-y_P0(:,i);
    deltU(:,i)=kmpc*(R-y_P0(:,i));
    U(:,i+1)=deltU(:,i)+U(:,i);
    y_N(:,i+1)=y_N0(:,i)+A_N*deltU(:,i);
    for j=1:p
        y(j,i+1)=y_N((j-1)*N+1,i+1);
    end
    t=0:delta:delta*i;
     y_Real1=lsim(sys(:,1),U(1,1:i+1),t)+lsim(sys(:,2),U(2,1:i+1),t);
     y_Real=y_Real1';
end

n=size(y_Real(1,:));
n=n(1,2);
t=1:1:n;

figure(2);
subplot(2,1,1);
plot(t,y_Real(1,1:n));
title('实际输出');
ylabel('y1','rotation',0);xlabel('time');

subplot(2,1,2);
plot(t,y_Real(2,1:n));
ylabel('y2','rotation',0);xlabel('time');

