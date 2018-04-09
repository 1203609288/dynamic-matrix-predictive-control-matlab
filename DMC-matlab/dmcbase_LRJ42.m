% 本文件演示最最基本的DMC算法
% 当 计算tk+1时,  滚动 y1(k+2/k) 生成 y0(k+1+1/k+1), 则y0(k+1+HP/k+1) ~= y1(k+HP/k) 会产生 y(HP)=y(HP)的问题, 为修正此问题, LRJ演示了四种修正方法

% LRJ修正方法1: 认为 阶跃响应是连续的平滑的,  可以认为是成功的, 使得 系统平稳 HP 个周期后的波动减小
%  DMC_y0 (m_nHP) += DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2) % 

% LRJ修正方法2:  加入死区，理论上可取，成功. 但是如果死区太大, 则势必使得系统的抗干扰性能降低
		%if( -0.6<duk & duk<0.6 )
		%	duk = 0 ;
		%end
		
% LRJ修正方法3: 不成功

% 	DMC_step(m_nHP) = y_dmc_step(m_nHP-1);		% 当计算 y1(k+HP/k) 时用的阶跃响应系数 修改了.
%   因为在计算 y1(k+1+HP/k+1) = y0(k+1+HP/k+1) + a(HP) * du(k+1) 会有误差
%   故采用     y1(k+1+HP/k+1) = y0(k+1+HP/k+1) + a(HP-1) * du(k+1), 是考虑到 其实  y0(k+1+HP/k+1)在数值上 =  y0(k+HP/k+1)
%   这是因为滚动的原因, y0(k+1+HP/k+1) = y0(k+HP/k+1), 因此 y1(k+1+HP/k+1) 实际在数值上= y1(k+HP/k+1)
%
% LRJ修正方法4: 这是第四种修正方法,也是不成功的.
%	y0(HP+1) 的长度为 HP+1,
%   而使用DMC控制的事后, 时域为 HP.
%   仿真结果表明这只是把小波动延后了一个周期,并无实质性的改变. 相当与 HP增加了, 控制效果好一些.
%
% 方法4 有两种实现方法,
% 4.1 是把 y0 设为(m_nHP+1,1),
% 4.2 是引入变量 y0p1
% 本文件采用方法4.2
clear;
clc;
bShowMsg = false;
% 第一部分：被控对象的参数,用于建立被控对象的模型,以及DMC控制的仿真.
% 被控对象参数Gp
T0_gp = [0.01, 0.3455,0.084, 8];
k0_gp = [1, 1, 1,20];
n_gp=4;
h_gp  = 0.5;
a_gp = exp( -h_gp ./ T0_gp );
b_gp = k0_gp .* ( 1 - a_gp );

Ts = 2;						% DMC闭环控制,采样周期
m_nHP = 10;					% DMC 参数, 预测时域
m_nHM = 3;					% DMC 参数, 控制时域

[ADMC, A, y_dmc_step, m_nHP, m_nHM ] = CalcStepResponseADMC( T0_gp, k0_gp, Ts, m_nHP, m_nHM , false);

m_nHPideal = length(y_dmc_step) ;	% 被控Gp达到稳态, 完整的预测时域,

DMC_step = y_dmc_step(1:m_nHP);
% DMC_step(m_nHP) = y_dmc_step(m_nHP-1);	% LRJ修正方法3: 不成功


%load dmc1_log ADMC y_dmc_step
%load dmc1_d0_log
%dModel=d_dmc_step(1:m_nHP);

% 第二部分：扰动通道
T0_gd=8;	k0_gd = 1;

h_gd =0.25;			% 扰动通道参数. 仿真时间步长=0.25
a_gd = exp(-h_gd / T0_gd); b_gd = k0_gd*(1-a_gd);


% 第三部分：被控对象的DMC控制, 

% 被控对象的输入输出
ysp = 1;	% 闭环控制, 设定输入信号
duk = 0;	% u(k) - u(k-1)
uk = 0 ;    % mv
yk = 0      % 被控对象的输出

dk=0.2;		% 扰动信号，是扰动环节的输入端
ydk=dk;		% 扰动环节的输出

itEnd	= 500;					% 仿真长度
t_time	= zeros(itEnd,1);		% 仿真时间坐标, 从0开始. t_time(1)=0, t_time(2)=h_gp, t_time(itEnd) = (itEnd-1)*h_gp

u_mv	= zeros(itEnd,1);		% uk, 存储历史数据,用于作图.
y_pv	= zeros(itEnd,1);		% yk	

d_dv	= zeros(itEnd,1);		% dk
y_dv	= zeros(itEnd,1);		% ydk

y_sum	= zeros(itEnd,1);		% yk + ydk, 

DMC_yideal	= zeros( size(y_dmc_step) );	% 到达阶跃稳态的响应, HP足够长, 理想状态

DMC_y0	= zeros( m_nHP, 1 );			% 无 deltaU(k) 作用时的自由响应, y0(k+i/k).		计算出deltaU(k)后,则被更新,相当于 y1(k+i/k) = y0(k+i/k+1)
DMC_ye	= zeros( m_nHP, 1 );			% 偏差, w(k+i/k) - y(k+i/k)
DMC_du	= zeros( m_nHM, 1 );			% deltaU(k), mv增量
DMC_yw	= ysp * ones( m_nHP, 1 );		% 参考轨迹. w(k+i/k),	本仿真简单处理=1

DMC_du_mem = zeros(itEnd, m_nHM);		%DMC计算过程的存储,用于分析计算过程
DMC_y0_mem = zeros(itEnd, m_nHP);
DMC_y01_mem = zeros(itEnd, m_nHP);		%计算出deltaU(k)后,则被更新,相当于 y1(k+i/k) = y0(k+i/k+1)
DMC_ye_mem = zeros(itEnd, m_nHP);		

DMC_yideal_mem = zeros(itEnd, m_nHPideal);

it = 1; itdmc=0;				% 循环变量
tk = 0; 
yerror = 0;						% 闭环校正, 实测误差 = y(k+1) - y1(k+1/k)

clc;
while (it < itEnd)
	it = it + 1;				% 第一次, it=2, y(2)=t1时刻, y(1)=0时刻的初值.
	tk = tk + h_gp;				% 当前时刻, 用于调试
	
	%----	 扰动环节的计算		----
	
	%	 扰动输入信号
	if( 50<it && it<70)     dk = 0;		% 如果不施加扰动，则令 dk = 0 
    else                    dk = 0.2;
	end
	
	% 扰动的步长 h_gd = 0.25, 而 被控对象的仿真步长 h_gp = 0.5, 因此需要计算两次
	ydk = a_gd * ydk + b_gp * dk;
	ydk = a_gd * ydk + b_gp * dk;
    
	ydk = 0.2;					% 如果不考虑扰动信号,则令 =0
	
	%---- 被控对象的仿真计算	----
	yk = a_gp * yk + b_gp * uk ;

    %if( yk <= 0.9 & it > 10 )
     %   it, uk, yk;
    %end
	
	yksum = yk + ydk;			% yksum = 被控对象的输出 y_pv + 扰动环节的输出 ydk
	
	if( mod(it, 4) == 1 )	
        
        itdmc = itdmc + 1;		% 第一次, itdmc=1
        
		if( itdmc == 1 )		% 第一次执行DMC控制
			DMC_y0 = yk* ones( size(DMC_y0) );	y0p1 = yk;
			DMC_yideal = yk* ones( size(DMC_yideal) );
		else
			yerror = yksum - DMC_y0(1) ;								% 实测偏差
			DMC_y0 ( 2 : m_nHP )   = DMC_y0 ( 2 : m_nHP ) + yerror; 	% 反馈校正	
			y0p1 = y0p1 + yerror;			
			
			DMC_y0 ( 1 : m_nHP-1 ) = DMC_y0 ( 2 : m_nHP );				% 滚动
			DMC_y0 ( m_nHP )	   = y0p1;
				
			% LRJ修正方法1: 认为 阶跃响应继续, 可以认为是成功的
			% DMC_y0 (m_nHP) = DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2)+DMC_y0(m_nHP); % 
			                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          % DMC_y0 (m_nHP) = DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2) + DMC_y0(m_nHP); 
			yerror 		= yksum - DMC_yideal(1);
			DMC_yideal ( 2 : m_nHPideal   ) = DMC_yideal ( 2 : m_nHPideal ) + yerror; 
			DMC_yideal ( 1 : m_nHPideal-1 ) = DMC_yideal ( 2 : m_nHPideal );
		end
		
		% DMC_yd = dModel * dk;					% 对扰动ydk的一步预测
		% 设定值 - 预测y  -前馈补偿
		% DMC_ye = DMC_yw - DMC_y0 - DMC_yd;	% 设定值-预测y-前馈补偿

		% 下式没有前馈补偿
		DMC_ye = DMC_yw - DMC_y0 ;		% 设定值-预测y		
		
		DMC_du = ADMC * DMC_ye ;		% 计算控制规律 Δuk		
		duk = DMC_du(1);
		
		% LRJ修正方法2:  加入死区，理论上可取，成功
		%if( -0.6<duk & duk<0.6 )
		%	duk = 0 ;
		%end
		uk = uk + duk;
		
		% PID 控制,  
		% 当K>2时, 会出现数值计算不稳定的问题,导致仿真结果震荡,发散.
		% uk = 10 * (ysp - yksum);
		
		%----  调试信息	---------------------------
		%----  调试信息,发现偏差
		if ( bShowMsg )
			msg = sprintf('it=%d, itdmc= %d, time=%g, \t yerror=%f, ' , it, itdmc, (it-1)*h_gp, yerror) ;
			disp(msg);

			% [[1:50]',DMC_y0, DMC_ye];		
			if( DMC_du(1) > 0.01 && it > 10)
				[u_mv(it-10:it), y_pv(it-10:it)]
				DMC_ye,
			end
		end
		%----  调试信息,发现偏差
          
        % DMC计算过程的存储,用于分析计算过程
     	DMC_ye_mem(itdmc,:) = DMC_ye';
        DMC_y0_mem(itdmc,:) = DMC_y0';	
		DMC_yideal_mem(itdmc,:) = DMC_yideal';
			
		DMC_du_mem(itdmc,:)  = DMC_du';   
		%----  调试信息	---------------------------
		
		 
		% DMC, delU(k) 起作用后，y1(k+1/k)
        DMC_y0 		= DMC_y0 	 + duk * DMC_step; 			% 计算出 y1
		y0p1		= y0p1       + duk * y_dmc_step(m_nHP+1);
		
		DMC_yideal	= DMC_yideal + duk * y_dmc_step;		% 计算出 y1,HP足够长
		
		%----  调试信息	---------------------------
		DMC_y01_mem(itdmc,:) = DMC_y0';			% 存储DMC计算过程
		%----  调试信息	---------------------------
		
    end
        % 保存到 数组中，用于后面的绘图
        d_dv(it)	= dk;		y_dv(it)	= ydk;		% 扰动通道 gd 的输入输出
        u_mv(it)	= uk;    	y_pv(it)	= yk;		% 被控对象 gp 的输入输出
        y_sum(it)	= yksum;

        t_time(it)	= tk;
end

figure();
msg = sprintf('没有前馈补偿的DMC Controler: HP=%d, HM=%d, A(%d,1)=%f.LRJ修正42', m_nHP, m_nHM, m_nHP,A(m_nHP,1));

% t_time = 0: h_gp : (itEnd-1)*h_gp;
	
subplot(211);
%plot(t_time, y_pv);
plot(t_time, y_pv, 'r-', t_time, y_dv, 'g-', t_time, y_sum, 'k-' );
ylabel('y');
legend('被控对象的输出ypv', '扰动输出yd', 'DMC的被调量(ypv+yd)');

title(msg);
grid on;

subplot(212);
stairs(t_time, u_mv, 'r-');	hold on;
% stairs(t_time, d_dv,'g-');
legend('DMC控制作用u_mv', '扰动输入 d_dv');
ylabel('u');
xlabel('time');
title('控制作用 u, 扰动作用 d, dead=0');
grid on;