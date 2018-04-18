function [D,A,L,A_N,a,Q] = DMC_martixD( model,uwt,ywt,M,P,m,p,N )
%本函数的目的是计算控制矩阵D

%计算R矩阵
R_t=[];
for j=1:1:m
    pj=uwt(1,j);
    Pj=pj*eye(M);
    R_t=blkdiag(R_t,Pj);
end

%   得到A矩阵
for i=1:1:p
    for j=1:1:m
        for l=1:1:N
            a(i,j,l)=model((l-1)*p+i,j);%输出yi对输入uj的响应
        end
    end
end

%得到Q矩阵 Q=eye(2*P);
Q=[];
for i=1:1:p
    qi=ywt(1,i);
    Qi=qi*eye(P);
    Q=blkdiag(Q,Qi);
end
% for i=1:1:3
%     QI((i-1)*15+15,(i-1)*15+15)=100;
% end;

for i=1:1:p
    for j=1:1:m
        for l=1:1:M
            A(((i-1)*P+l):(i*P),(j-1)*M+l)=a(i,j,1:(P-l+1));
        end
    end
end

for i=1:1:p
    for j=1:1:m
        A_N((i-1)*N+1:i*N,j)=a(i,j,:);
    end
end


%   得到L矩阵
for i=1:1:m
    L(i,(i-1)*M+1)=1;
    for j=2:1:M
        L(i,(i-1)*M+j)=0;
    end
end

%计算控制器增益矩阵
D=L*inv(A'*Q*A+R_t)*A'*Q; %课本第103页(6-8)

end

