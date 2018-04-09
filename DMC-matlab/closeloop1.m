
%  闭环控制， 控制器 纯比例控制器， Kp =1 
%
%

T0_gp=20  ;	
k0_gp =10;
Ts = 2;

h_gp = T0_gp / 50;			% 被控对象 仿真 时间步长
a_gp = exp(-h_gp / T0_gp); b_gp = k0_gp*(1-a_gp);

uk = 0 ;
yk = 0;
ysp = 1;

dk = 0;

it = 0;
itEnd = 200;
while (it < itEnd)
	it = it + 1;
	
	if( 50<it && it<70)     dk = 5;		% 如果不施加扰动，则令 dk = 0 
    else                    dk = 0;
	end

    ydk = 0; % dk;
	
	
	% 被控对象仿真
	yk = a_gp * yk + b_gp * uk ;
	
	yksum = yk + ydk;
    
    
    error = ysp - yksum;
    uk = error;
    
    y_d(it) = ydk;  
    
    u_mv(it) = uk;
    y_pv(it) = yk;    
    y_sum(it) = yksum;		% y_sum = 被控对象的输出 y_pv + 扰动环节的输出 ydk
    
end

figure();
t_time = h_gp : h_gp :  itEnd*h_gp;

subplot(211);
plot( t_time, y_pv, 'r', t_time, y_sum, 'b') ;
subplot(212);
plot(t_time, u_mv, 'r');



    