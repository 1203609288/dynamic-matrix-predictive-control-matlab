clear;clc;
num={[5],[3];[6],[9]};
den={[3 1 3],[1 2 5];[2 1 7],[2 3 6]};
sys=tf(num,den);%ģ�ʹ��ݺ���
g11=poly2tfd(num{1,1},den{1,1},0,0);
g12=poly2tfd(num{1,2},den{1,2},0,0);
g21=poly2tfd(num{2,1},den{2,1},0,0);
g22=poly2tfd(num{2,2},den{2,2},0,0);
delta=0.5;%����ʱ��
P=12;M=6;m=2;p=2;N=40; %   M,P,m,p�ֱ�Ϊ����ʱ�򳤶ȣ�Ԥ��ʱ�򳤶ȣ������������������������NΪ��ģʱ��
ny=2;
tfinal=500;
mymodel=tfd2step(tfinal,delta,ny,g11,g12,g21,g22);%�����Ծ��Ӧģ��

%��ͼ������ģ�͵Ľ�Ծ��Ӧ����
figure(1)
subplot(2,2,1);
step(num{1,1},den{1,1});
title('u1-y1��Ծ��Ӧ');
xlabel('time');
subplot(2,2,2);
step(num{1,2},den{1,2});
title('u1-y2��Ծ��Ӧ');
xlabel('time');
subplot(2,2,3);
step(num{2,1},den{2,1});
title('u2-y1��Ծ��Ӧ');
xlabel('time');
subplot(2,2,4);
step(num{2,2},den{2,2});
title('u2-y2��Ծ��Ӧ');
xlabel('time');

ywt=[3,1];%Q����
uwt=[400,300];%R����
alpha=[1,1];%H����
r=[1;2];%�趨ֵ
tend=500;%����ʱ��

%���㷴��У��H����
H=[];
for i=1:p
    h=alpha(1,i)*ones(N,1);
    H=blkdiag(H,h);
end

%������λ����S
for i=1:p
    for j=1:N-1
        S((i-1)*N+j,(i-1)*N+j+1)=1;%�ζԽ���Ԫ��Ϊ1
    end
    S((i*N),(i*N))=1;%���½�Ԫ��Ϊ1
end

%����趨ֵR����
R=[];
for i=1:p
    r=r(i,1)*ones(P,1);
    R=[R;r];
end

y_Real=zeros(p,tend);%ʵ�����
e=zeros(p,tend);%���
y=zeros(p,tend);%����
U=zeros(m,tend);%����

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


[kmpc,A,L,A_N,a,Q]=DMC_martixD(mymodel,uwt,ywt,M,P,m,p,N);%����D����

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
title('ʵ�����');
ylabel('y1','rotation',0);xlabel('time');

subplot(2,1,2);
plot(t,y_Real(2,1:n));
ylabel('y2','rotation',0);xlabel('time');

