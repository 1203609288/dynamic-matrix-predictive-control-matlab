clear;
clc;
% 第一部分：被控对象的仿真
% 被控对象参数
T0_gp=25;	k0_gp = 10;

h_gp = T0_gp / 50;			% 被控对象 仿真 时间步长 0.5
a_gp = exp(-h_gp / T0_gp); b_gp = k0_gp*(1-a_gp);

Ts = 2;
m_nHP = 10;
m_nHM = 3;

[ADMC, A, y_dmc_step, m_nHP, m_nHM ] = CalcStepResponseADMC( T0_gp, k0_gp, Ts, m_nHP, m_nHM , false);
m_nHP0 = length(y_dmc_step) ;

%load dmc1_log ADMC y_dmc_step
%load dmc1_d0_log
%dModel=d_dmc_step(1:m_nHP);

% 第二部分：扰动通道
T0_gd=10;	k0_gd = 5;

h_gd = T0_gd / 40;			% 扰动通道 仿真 时间步长， 0.25
a_gd = exp(-h_gd / T0_gd); b_gd = k0_gd*(1-a_gd);



% 第三部分：被控对象的DMC控制, 

% 被控对象的输入输出
ysp = 1;	% 设定输入信号
uk = 0 ;    % mv
yk = 0      % 被控对象的输入和输出

dk=0;		% 扰动信号，是扰动环节的输入端
ydk=0;		% 扰动环节的输出

itEnd = 500;
t_time =  zeros(itEnd,1);

u_mv =  zeros(itEnd,1);			% 存储历史数据，用于作图
y_pv =  zeros(itEnd,1);			
y_d  = zeros(itEnd,1);
y_sum = zeros(itEnd,1);

DMC_y1 = zeros( size(y_dmc_step) );	% 到达阶跃稳态的响应,HP足够长

DMC_y0 = zeros( m_nHP, 1 );		% 无 deltalU(k) 作用时的自由响应
DMC_ye = zeros( m_nHP, 1 );		% 偏差
DMC_du = zeros( m_nHM, 1 );		% mv增量
DMC_yw = ysp * ones( m_nHP, 1 );		% 参考轨迹. 本仿真简单处理=1

DMC_du_mem = zeros(itEnd, m_nHM);	%DMC计算过程的存储
DMC_y0_mem = zeros(itEnd, m_nHP);
DMC_y1_mem = zeros(itEnd, m_nHP0);

it = 1; itdmc=0;
yerror = 0;

clc;
while (it < itEnd)
	it = it + 1;
	
	%----	 扰动环节的计算		----
	
	%	 扰动输入信号
	if( 50<it && it<70)     dk = 5;		% 如果不施加扰动，则令 dk = 0 
    else                    dk = 0;
	end
	
	% 扰动的步长 h_gd = 0.25, 而 被控对象的仿真步长 h_gp = 0.5 
	% 因此需要计算2次
	ydk = a_gd * ydk + b_gp * dk;
	ydk = a_gd * ydk + b_gp * dk;
    
	ydk = 0;
	
	%---- 被控对象的仿真计算	----
	yk = a_gp * yk + b_gp * uk ;

    if( yk <= 0.9 & it > 10 )
        it, uk, yk;
    end
	
	% y_sum = 被控对象的输出 y_pv + 扰动环节的输出 ydk
	yksum = yk + ydk;
	
	if( mod(it, 4) == 1 )	
        
        itdmc = itdmc + 1;
        
		if( itdmc == 1 )	% 第一次执行DMC控制
			DMC_y0 = yk*ones(size(DMC_y0));
			DMC_y1 = yk*ones(size(DMC_y1));
		else
			yerror = yksum - DMC_y0(1) ;			% 偏差
			DMC_y0 ( 2 : m_nHP )   = DMC_y0 ( 2 : m_nHP ) + yerror; 	% 反馈校正				
			DMC_y0 ( 1 : m_nHP-1 ) = DMC_y0 ( 2 : m_nHP );			% 滚动
				
			% DMC_y0 (m_nHP) = DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2)+DMC_y0(m_nHP); % LRJ的矫正

			DMC_y1 ( 2 : m_nHP0 )   = DMC_y1 ( 2 : m_nHP0 ) + yerror; 
			DMC_y1 ( 1 : m_nHP0-1 ) = DMC_y1 ( 2 : m_nHP0 );
		end
		
		% DMC_yd = dModel * dk;				% 对扰动ydk的一步预测
		
		% 设定值-预测y-前馈补偿
		% DMC_ye = DMC_yw - DMC_y0 - DMC_yd;	% 设定值-预测y-前馈补偿

		% 下式没有前馈补偿
		DMC_ye = DMC_yw - DMC_y0 ;		% 设定值-预测y		
		
		DMC_du = ADMC * DMC_ye ;		% 计算控制规律 Δuk			
		uk = uk + DMC_du(1);
		
		% PID 控制,  当K>2时, 会出现数值计算不稳定的问题,导致仿真结果震荡,发散.
		% uk = 10 * (ysp - yksum);  
        msg = sprintf('it=%d, itdmc= %d , time =%g, yerror=%f, ' , it, itdmc, (it-1)*h_gp, yerror) ;
		disp(msg);
          
        DMC_du_mem(itdmc,:)  = DMC_du';
        DMC_y0_mem(itdmc,:) = DMC_y0';
		DMC_y1_mem(itdmc,:) = DMC_y1';
        
        % [[1:50]',DMC_y0, DMC_ye];		
        if( DMC_du(1) > 0.01 && it > 10)
            [u_mv(it-10:it), y_pv(it-10:it)]
			DMC_ye,
		end
		 
        DMC_y0 = DMC_y0 + DMC_du(1) * A(1:m_nHP,1); 	% 计算出 y1
		DMC_y1 = DMC_y1 + DMC_du(1) * y_dmc_step;		% 计算出 y1,HP足够长      
		
    end
        % 保存到 数组中，用于 后面的绘图
        u_mv(it) = uk;    
        y_d (it) = ydk;
        y_pv(it) = yk;	
        y_sum(it)= yksum;    
end

figure();
msg = sprintf('没有前馈补偿的DMC Controler: HP=%d, HM=%d, A(%d,1)=%f', m_nHP, m_nHM, m_nHP,A(m_nHP,1));

subplot(211);
t_time = 0: h_gp : (itEnd-1)*h_gp;
%plot(t_time, y_pv);
plot(t_time, y_pv, 'r-', t_time, y_d, 'g-', t_time, y_sum, 'k-' );
ylabel('y');
title(msg);
legend('被控对象的输出ypv', '扰动输出ykd', 'DMC的被调量(ypv+ykd)');
grid on;

subplot(212);
stairs(t_time, u_mv);
ylabel('u');
xlabel('time');
title('控制作用 uk');
grid on;
