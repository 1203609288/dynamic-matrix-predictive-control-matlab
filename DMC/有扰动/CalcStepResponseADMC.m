function [ADMC, A, y_step, m_nHP, m_nHM ] =  CalcStepResponseADMC( T0_gp, k0_gp, Ts, m_nHP, m_nHM , bFigure)

% 被控对象 仿真，获取 阶跃响应 向量
% y_dmc_stepy = 被控对象 到 稳态 的完整的阶跃响应
% A		: 动态矩阵
% ADMC	: deltaU = ADMC * DMC_ye
% m_nHP	: 预测时域长度
% m_nHM	: 控制时域长度

% 第一部分：被控对象的仿真
% 被控对象参数
% 
% 传递来的缺省参数
%T0_gp=25;	k0_gp = 10;

h_gp = T0_gp / 50;			% 被控对象 仿真 时间步长
a_gp = exp(-h_gp / T0_gp); b_gp = k0_gp*(1-a_gp);
Tend =  10 * T0_gp / h_gp;

y_pv = zeros(Tend, 1);			% 存储仿真结果,以h为仿真步长
y_step = zeros(m_nHP, 1);	% 控制器用,以 Ts 为采样周期

uk = 1;
y_pv(1) = 0;                % 初始状态 y(0)=0
yss = 0.9999*k0_gp * uk;

for ( it = 1 : Tend )
	y_pv(it+1) = a_gp * y_pv(it) + b_gp * uk;
	if( y_pv(it+1) > yss )
        itEnd = it+1;
        break;
    end
end

% size(t_time), size(y_pv),
if( bFigure )
    t_time = 0: h_gp: (itEnd-1)*h_gp;
    figure();
    plot(t_time, y_pv(1:itEnd), 'r-');  
    hold on;
end


% 第二部分：根据控制器的采样周期，生成 控制器用的 阶跃响应

% Ts = 2;			% 采样周期 =2 
mulTime = Ts / h_gp;

it2 = 1;			% 阶跃响应的 循环 a(it2)
it  = 1+mulTime;	% 4 = Ts / h 

while( it <= itEnd )
    y_step(it2) = y_pv(it);
    it = it + mulTime;
    it2 = it2 + 1;
end
% y_step ( it2 ) = y_pv(itEnd);     /// 是否 aN = ass

% 第三部分：计算DMC所用的 模型矩阵 A 
A=zeros(m_nHP, m_nHM);	% DMC的模型矩阵 A 
for( it = 1 : m_nHM )
	A([it : m_nHP], it) = y_step( [1 : m_nHP+1-it] );
end

if( bFigure )
    t2_time = Ts : Ts : (it2-1)*Ts;
	plot(t2_time, y_step, 'b.');
    plot(t2_time([1:m_nHP]), A([1:m_nHP],1), 'm*');
    grid on;
    legend ( 'y.pv(h)', 'y.step(Ts)', 'y.dmc.step(m.nHP*Ts)' );

	xlabel('time');
	ylabel('y');

	msg = sprintf('被控对象=%g/(%gs+1), h=%g, m_nHP=%d, m_nHM=%d', k0_gp, T0_gp, h_gp, m_nHP, m_nHM );
	title(msg,'Interpreter', 'none');
end

%if( m_nHP < length(y_step ) )
%    A(m_nHP, 1) = 10; % //y_step( length(y_step) );
%end

% 第四部分：计算DMC所用的 模型矩阵 A_DMC 

% DMC 的可整定的控制参数
m_fQ = eye( m_nHP, m_nHP );
% for( i=1:m_nHP )     m_fQ(i,i)= i*i; end
m_fR = eye( m_nHM, m_nHM ) * 0.2;
m_fR = zeros( m_nHM, m_nHM ) ;

% 下面，计算 DMC 矩阵 DMC = inv(A'QA+R) A'Q
At = A';
At1 = At * m_fQ;
At2 = At1 * A;
At3 = At2 + m_fR;
At4 = inv(At3);
At5 = At4*At1;
ADMC = At5

%A(m_nHP, 1) = 10;
% save dmc1_log ADMC y_step

% y_step ( m_nHP+1 , 1 ) = uk * k0_gp ;                % = 10, 研究ass的作用

save CalcStepResponseADMC y_step A ADMC