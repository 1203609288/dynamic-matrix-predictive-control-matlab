function [ADMC, A, y_step, m_nHP, m_nHM ] =  CalcStepResponseADMC( T0_gp, k0_gp, Ts, m_nHP, m_nHM , bFigure, nLoop, Kp2)

% 被控对象 仿真，获取 阶跃响应 向量
% y_step = 被控对象 到 稳态 的完整的阶跃响应
% A		: 动态矩阵
% ADMC	: deltaU = ADMC * DMC_ye
% m_nHP	: 预测时域长度
% m_nHM	: 控制时域长度

%m_nHP = 100;
%m_nHM = 5;
%bFigure = true;

% 第一部分：被控对象的仿真
% 被控对象参数
% 
% 传递来的缺省参数
% T0_gp=25;	k0_gp = 10;
% h_gp = T0_gp / 50;			% 被控对象 仿真 时间步长
% Tend =  10 * T0_gp / h_gp;

%---- 2013 试验汽温对象
%nLoop = 1;		% 内回路
%nLoop = 2;		% 外回路

n_gp=5;							% 5阶对象, 2阶在内回路，3阶在外回路
%T0_gp = [15, 15,	20,     20, 20];
%k0_gp = [8,  1,		1.125,  1,  1 ];
h_gp  = 0.5;					% 被控对象 仿真 时间步长

a_gp = exp( - h_gp ./ T0_gp );
b_gp = k0_gp .* ( 1 - a_gp );
iTimeSimuEnd = 4500;			% 仿真长度,
iTimeSS = 0;					% 达到稳态, static state

% y_pv(:, [内回路导前区汽温， 外回路惰性区汽温] )
y_pv = zeros(iTimeSimuEnd, 2);	% 存储仿真结果,以h为仿真步长
y_step = zeros(m_nHP,1);	    % 控制器用,以 Ts 为采样周期
y_pv(1,1) = 0;                % 初始状态 y1(0)=0
y_pv(1,2) = 0;                % 初始状态 y2(0)=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%导入采样数据值%%%%%%%%%%%%%%%%%%%%%%%%%%%

uk = 1;
if ( nLoop ==1 )		% 只内回路 的
	yss = 0.9999 * k0_gp(1) *k0_gp(2) * uk;		% 导前汽温稳态值
	nypvCol = 1;
elseif( nLoop == 2 )	% 只外回路 的
	yss = 0.9999 * k0_gp(1) *k0_gp(2) * k0_gp(3) *k0_gp(4) * k0_gp(5) * uk;	% 出口惰性区汽温
	nypvCol = 2;
elseif( nLoop == 3 )	% 内回路 P， 外回路 DMC
	Kinss = ( k0_gp(1) * k0_gp(2) * Kp2 );	
	Kinss = Kinss / (1 + Kinss );
	% yss = 0.9999 * Kinss / (1 + Kinss ) * uk ;	 % 以 内回路的稳态值，作为 ss
	yss = 0.9999 * Kinss *  k0_gp(3) *k0_gp(4) * k0_gp(5) * uk;	
	nypvCol = 2;
end	

x = zeros(n_gp+1,1);		% 仿真中间变量 x
 x(1) = uk;
for ( it = 1 : iTimeSimuEnd-1 )	
	if( nLoop ==1 | nLoop== 2)			
		x(1) = uk;			% 内回路开环的情形
	elseif( nLoop == 3 )	
		x(1) = Kp2 *( uk - x(3) ) ;		% 内回路采用P控制，KP=Kp2
	end
	
	for( inh = 1 : 5 )
		x(inh+1) = a_gp(inh) * x(inh+1) + b_gp(inh) * x(inh);
	end	
	y_pv(it+1,1) = x(3);	y_pv(it+1,2) = x(6);
	
%	达到稳态，则退出仿真
%	if( y_pv(it+1,1) > yss2 )
%       iTimeSS = it+1;
%        break;
%	end	
	if( (y_pv(it+1,nypvCol) > yss)  & (iTimeSS == 0 ) )
        iTimeSS = it+1;
        break;
	end	
end
y_len = max ( it+1, iTimeSS) ;	% 画图时，尽量取更多的数据

if( bFigure )
    t_time = 0: h_gp: (y_len-1)*h_gp;
    figure();
    plot(t_time, y_pv(1:y_len,1), 'r-', t_time, y_pv(1:y_len, 2), 'k-');
	% legend('导前区汽温，内回路', '惰性区汽温，外回路', 'Location', 'SouthEast');	
	grid on;
    hold on;
end


% 第二部分：根据控制器的采样周期，生成 控制器用的 阶跃响应
% Ts = 2;			% 采样周期 =2 
mulTime = Ts / h_gp;

it2 = 1;			% 阶跃响应的 循环 a(it2)
it  = 1+mulTime;	% 4 = Ts / h 
while( it <= iTimeSS )
    y_step(it2) = y_pv(it,nypvCol);		
	
    it = it + mulTime;
    it2 = it2 + 1;
end
% y_step ( it2 ) = y_pv(iTimeSS);     /// 是否 aN = ass

m_nHP = min(m_nHP, it2-1);

% 第三部分：计算DMC所用的 模型矩阵 A 
A=zeros(m_nHP, m_nHM);	% DMC的模型矩阵 A 
for( it = 1 : m_nHM )
	A([it : m_nHP], it) = y_step( [1 : m_nHP+1-it] );
end

if( bFigure)
    t2_time = Ts : Ts : (it2-1)*Ts;
	plot(t2_time, y_step(1 : length(t2_time) ), 'b.');
    plot(t2_time([1:m_nHP]), A([1:m_nHP],1), 'm*');
    grid on;
	legend('导前区汽温，内回路', '惰性区汽温，外回路', 'y.step(Ts)', 'y.dmc.step(m.nHP*Ts)', 'Location', 'SouthEast');
	
	xlabel('time');
	ylabel('y');
% 	msg = sprintf('被控对象: nLoop==%d, h=%g, Ts=%g, m_nHP=%d, m_nHM=%d',nLoop, h_gp, Ts, m_nHP, m_nHM );
% 	title(msg,'Interpreter', 'none');
title('再热部分阶跃响应变化');
end

%if( m_nHP < length(y_step ) )
%    A(m_nHP, 1) = 10; % //y_step( length(y_step) );
%end

% 第四部分：计算DMC所用的 模型矩阵 A_DMC 

% DMC 的可整定的控制参数
m_fQ = eye( m_nHP, m_nHP );
% for( i=1:m_nHP )     m_fQ(i,i)= i*i; end
m_fR = 5*eye( m_nHM, m_nHM );
% m_fR = zeros( m_nHM, m_nHM ) ;

% 下面，计算 DMC 矩阵 DMC = inv(A'QA+R) A'Q
At = A';				% A 的转置'
At1 = At * m_fQ;
At2 = At1 * A;
At3 = At2 + m_fR;
At4 = inv(At3);
At5 = At4*At1;
ADMC = At5

%A(m_nHP, 1) = 10;
% save dmc1_log ADMC y_step

% y_step ( m_nHP+1 , 1 ) = uk * k0_gp ;                % = 10, 研究ass的作用

if( nLoop==1)
	save CalcStepResponseADMC_in  y_step A ADMC
else
	save CalcStepResponseADMC_out y_step A ADMC
end
	