% 扰动通道 仿真，获取 阶跃响应 向量
% y_step     = 被控对象 到 稳态 的完整的阶跃响应
% y_dmc_step = y_dmc_stepy( 1 : m_nHP )
function [ y_dmc_step, y_step ] =  CalcStepResponse_D( T0_gd, k0_gd, Ts, m_nHP, m_nHM , bFigure)


disp('扰动通道的仿真');

% 第一部分：扰动通道的仿真

% 扰动通道参数
% T0_gd=10;		k0_gd = 5;
h_gd = T0_gd / 50;			% 仿真 时间步长,  0.2
a_gd = exp(-h_gd / T0_gd); b_gd = k0_gd*(1-a_gd);
iTend = 10 * T0_gd / h_gd;

y_pv		= zeros(iTend, 1);		% 存储仿真结果,以h为仿真步长
t_time		= zeros(iTend, 1);

y_step		= zeros(m_nHP, 1);	% 控制器用,以 Ts 为采样周期
y_dmc_step	= zeros(m_nHP, 1);	% 控制器用,以 Ts 为采样周期
uk = 1;

for ( it = 1 : iTend-1 )	% while ( it < Tend )
	y_pv(it+1) = a_gd * y_pv(it) + b_gd * uk;	% y_pv(1) = t0, y_pv(2) = t1

	% t_time(it+1) = t_time(it) + h_gd;			% 数值计算的误差更大一些
	t_time(it+1) = it * h_gd;

	%if( y_pv(it+1) > 0.98*k0_gd )
     %   itEnd = it+1;
     %  break;
    % end
end

itEnd = it+1;
t_time2 = [ 0: h_gd: (itEnd-1)*h_gd ]';				% 直接生成时间坐标，可能速度更快
% assert( t_time == t_time2 );						% 数值计算的误差，导致不相等，
if( bFigure)		% 如果需要，作图绘制 扰动的阶跃响应 （完整的）
	figure();
	plot(t_time, y_pv(1:itEnd), 'r-'); 		hold on;
	end

% 第二部分 根据控制器的采样周期，生成 控制器用的 阶跃响应
mulTime = Ts / h_gd;
it2 = 1;
it  = 1 + mulTime;

t_time2 = zeros(m_nHP, 1);
while( it <= itEnd )
    y_step(it2) = y_pv(it);
    t_time2(it2) = t_time(it);
    it = it + mulTime;
    it2 = it2 + 1;
end

y_dmc_step = y_step(1 : m_nHP);			% 只要求的控制时域的坐标
% y_dmc_step ( it2 ) = y_pv(itEnd);

% t_time3 = Ts : Ts : (it2-1)*Ts;
% assert( t_time3 == t_time2 );

if( bFigure )
	plot(t_time2, y_step, 'b.');
	plot(t_time2(1:m_nHP), y_dmc_step, 'm*');
	grid on;

	legend ( 'y.pv(h)', 'y.step(Ts)', 'y.dmc.step(m_nHP*Ts)' );

	xlabel('time');
	ylabel('y');

	msg = sprintf('扰动对象=%g/(%gs+1), h=%g, m_nHP=%d', k0_gd, T0_gd, h_gd, m_nHP );
	title(msg,'Interpreter', 'none');


end

save CalcStepResponse_D y_dmc_step y_step