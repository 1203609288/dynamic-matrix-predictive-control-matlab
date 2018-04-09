% ���ļ���ʾ���������DMC�㷨
% �� ����tk+1ʱ,  ���� y1(k+2/k) ���� y0(k+1+1/k+1), ��y0(k+1+HP/k+1) ~= y1(k+HP/k) ����� y(HP)=y(HP)������, Ϊ����������, LRJ��ʾ��������������

% LRJ��������1: ��Ϊ ��Ծ��Ӧ��������ƽ����,  ������Ϊ�ǳɹ���, ʹ�� ϵͳƽ�� HP �����ں�Ĳ�����С
%  DMC_y0 (m_nHP) += DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2) % 

% LRJ��������2:  ���������������Ͽ�ȡ���ɹ�. �����������̫��, ���Ʊ�ʹ��ϵͳ�Ŀ��������ܽ���
		%if( -0.6<duk & duk<0.6 )
		%	duk = 0 ;
		%end
		
% LRJ��������3: ���ɹ�

% 	DMC_step(m_nHP) = y_dmc_step(m_nHP-1);		% ������ y1(k+HP/k) ʱ�õĽ�Ծ��Ӧϵ�� �޸���.
%   ��Ϊ�ڼ��� y1(k+1+HP/k+1) = y0(k+1+HP/k+1) + a(HP) * du(k+1) �������
%   �ʲ���     y1(k+1+HP/k+1) = y0(k+1+HP/k+1) + a(HP-1) * du(k+1), �ǿ��ǵ� ��ʵ  y0(k+1+HP/k+1)����ֵ�� =  y0(k+HP/k+1)
%   ������Ϊ������ԭ��, y0(k+1+HP/k+1) = y0(k+HP/k+1), ��� y1(k+1+HP/k+1) ʵ������ֵ��= y1(k+HP/k+1)
%
% LRJ��������4: ���ǵ�������������,Ҳ�ǲ��ɹ���.
%	y0(HP+1) �ĳ���Ϊ HP+1,
%   ��ʹ��DMC���Ƶ��º�, ʱ��Ϊ HP.
%   ������������ֻ�ǰ�С�����Ӻ���һ������,����ʵ���Եĸı�. �൱�� HP������, ����Ч����һЩ.
%
% ����4 ������ʵ�ַ���,
% 4.1 �ǰ� y0 ��Ϊ(m_nHP+1,1),
% 4.2 ��������� y0p1
% ���ļ����÷���4.2
clear;
clc;
bShowMsg = false;
% ��һ���֣����ض���Ĳ���,���ڽ������ض����ģ��,�Լ�DMC���Ƶķ���.
% ���ض������Gp
T0_gp = [0.01, 0.3455,0.084, 8];
k0_gp = [1, 1, 1,20];
n_gp=4;
h_gp  = 0.5;
a_gp = exp( -h_gp ./ T0_gp );
b_gp = k0_gp .* ( 1 - a_gp );

Ts = 2;						% DMC�ջ�����,��������
m_nHP = 10;					% DMC ����, Ԥ��ʱ��
m_nHM = 3;					% DMC ����, ����ʱ��

[ADMC, A, y_dmc_step, m_nHP, m_nHM ] = CalcStepResponseADMC( T0_gp, k0_gp, Ts, m_nHP, m_nHM , false);

m_nHPideal = length(y_dmc_step) ;	% ����Gp�ﵽ��̬, ������Ԥ��ʱ��,

DMC_step = y_dmc_step(1:m_nHP);
% DMC_step(m_nHP) = y_dmc_step(m_nHP-1);	% LRJ��������3: ���ɹ�


%load dmc1_log ADMC y_dmc_step
%load dmc1_d0_log
%dModel=d_dmc_step(1:m_nHP);

% �ڶ����֣��Ŷ�ͨ��
T0_gd=8;	k0_gd = 1;

h_gd =0.25;			% �Ŷ�ͨ������. ����ʱ�䲽��=0.25
a_gd = exp(-h_gd / T0_gd); b_gd = k0_gd*(1-a_gd);


% �������֣����ض����DMC����, 

% ���ض�����������
ysp = 1;	% �ջ�����, �趨�����ź�
duk = 0;	% u(k) - u(k-1)
uk = 0 ;    % mv
yk = 0      % ���ض�������

dk=0.2;		% �Ŷ��źţ����Ŷ����ڵ������
ydk=dk;		% �Ŷ����ڵ����

itEnd	= 500;					% ���泤��
t_time	= zeros(itEnd,1);		% ����ʱ������, ��0��ʼ. t_time(1)=0, t_time(2)=h_gp, t_time(itEnd) = (itEnd-1)*h_gp

u_mv	= zeros(itEnd,1);		% uk, �洢��ʷ����,������ͼ.
y_pv	= zeros(itEnd,1);		% yk	

d_dv	= zeros(itEnd,1);		% dk
y_dv	= zeros(itEnd,1);		% ydk

y_sum	= zeros(itEnd,1);		% yk + ydk, 

DMC_yideal	= zeros( size(y_dmc_step) );	% �����Ծ��̬����Ӧ, HP�㹻��, ����״̬

DMC_y0	= zeros( m_nHP, 1 );			% �� deltaU(k) ����ʱ��������Ӧ, y0(k+i/k).		�����deltaU(k)��,�򱻸���,�൱�� y1(k+i/k) = y0(k+i/k+1)
DMC_ye	= zeros( m_nHP, 1 );			% ƫ��, w(k+i/k) - y(k+i/k)
DMC_du	= zeros( m_nHM, 1 );			% deltaU(k), mv����
DMC_yw	= ysp * ones( m_nHP, 1 );		% �ο��켣. w(k+i/k),	������򵥴���=1

DMC_du_mem = zeros(itEnd, m_nHM);		%DMC������̵Ĵ洢,���ڷ����������
DMC_y0_mem = zeros(itEnd, m_nHP);
DMC_y01_mem = zeros(itEnd, m_nHP);		%�����deltaU(k)��,�򱻸���,�൱�� y1(k+i/k) = y0(k+i/k+1)
DMC_ye_mem = zeros(itEnd, m_nHP);		

DMC_yideal_mem = zeros(itEnd, m_nHPideal);

it = 1; itdmc=0;				% ѭ������
tk = 0; 
yerror = 0;						% �ջ�У��, ʵ����� = y(k+1) - y1(k+1/k)

clc;
while (it < itEnd)
	it = it + 1;				% ��һ��, it=2, y(2)=t1ʱ��, y(1)=0ʱ�̵ĳ�ֵ.
	tk = tk + h_gp;				% ��ǰʱ��, ���ڵ���
	
	%----	 �Ŷ����ڵļ���		----
	
	%	 �Ŷ������ź�
	if( 50<it && it<70)     dk = 0;		% �����ʩ���Ŷ������� dk = 0 
    else                    dk = 0.2;
	end
	
	% �Ŷ��Ĳ��� h_gd = 0.25, �� ���ض���ķ��沽�� h_gp = 0.5, �����Ҫ��������
	ydk = a_gd * ydk + b_gp * dk;
	ydk = a_gd * ydk + b_gp * dk;
    
	ydk = 0.2;					% ����������Ŷ��ź�,���� =0
	
	%---- ���ض���ķ������	----
	yk = a_gp * yk + b_gp * uk ;

    %if( yk <= 0.9 & it > 10 )
     %   it, uk, yk;
    %end
	
	yksum = yk + ydk;			% yksum = ���ض������� y_pv + �Ŷ����ڵ���� ydk
	
	if( mod(it, 4) == 1 )	
        
        itdmc = itdmc + 1;		% ��һ��, itdmc=1
        
		if( itdmc == 1 )		% ��һ��ִ��DMC����
			DMC_y0 = yk* ones( size(DMC_y0) );	y0p1 = yk;
			DMC_yideal = yk* ones( size(DMC_yideal) );
		else
			yerror = yksum - DMC_y0(1) ;								% ʵ��ƫ��
			DMC_y0 ( 2 : m_nHP )   = DMC_y0 ( 2 : m_nHP ) + yerror; 	% ����У��	
			y0p1 = y0p1 + yerror;			
			
			DMC_y0 ( 1 : m_nHP-1 ) = DMC_y0 ( 2 : m_nHP );				% ����
			DMC_y0 ( m_nHP )	   = y0p1;
				
			% LRJ��������1: ��Ϊ ��Ծ��Ӧ����, ������Ϊ�ǳɹ���
			% DMC_y0 (m_nHP) = DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2)+DMC_y0(m_nHP); % 
			                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          % DMC_y0 (m_nHP) = DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2) + DMC_y0(m_nHP); 
			yerror 		= yksum - DMC_yideal(1);
			DMC_yideal ( 2 : m_nHPideal   ) = DMC_yideal ( 2 : m_nHPideal ) + yerror; 
			DMC_yideal ( 1 : m_nHPideal-1 ) = DMC_yideal ( 2 : m_nHPideal );
		end
		
		% DMC_yd = dModel * dk;					% ���Ŷ�ydk��һ��Ԥ��
		% �趨ֵ - Ԥ��y  -ǰ������
		% DMC_ye = DMC_yw - DMC_y0 - DMC_yd;	% �趨ֵ-Ԥ��y-ǰ������

		% ��ʽû��ǰ������
		DMC_ye = DMC_yw - DMC_y0 ;		% �趨ֵ-Ԥ��y		
		
		DMC_du = ADMC * DMC_ye ;		% ������ƹ��� ��uk		
		duk = DMC_du(1);
		
		% LRJ��������2:  ���������������Ͽ�ȡ���ɹ�
		%if( -0.6<duk & duk<0.6 )
		%	duk = 0 ;
		%end
		uk = uk + duk;
		
		% PID ����,  
		% ��K>2ʱ, �������ֵ���㲻�ȶ�������,���·�������,��ɢ.
		% uk = 10 * (ysp - yksum);
		
		%----  ������Ϣ	---------------------------
		%----  ������Ϣ,����ƫ��
		if ( bShowMsg )
			msg = sprintf('it=%d, itdmc= %d, time=%g, \t yerror=%f, ' , it, itdmc, (it-1)*h_gp, yerror) ;
			disp(msg);

			% [[1:50]',DMC_y0, DMC_ye];		
			if( DMC_du(1) > 0.01 && it > 10)
				[u_mv(it-10:it), y_pv(it-10:it)]
				DMC_ye,
			end
		end
		%----  ������Ϣ,����ƫ��
          
        % DMC������̵Ĵ洢,���ڷ����������
     	DMC_ye_mem(itdmc,:) = DMC_ye';
        DMC_y0_mem(itdmc,:) = DMC_y0';	
		DMC_yideal_mem(itdmc,:) = DMC_yideal';
			
		DMC_du_mem(itdmc,:)  = DMC_du';   
		%----  ������Ϣ	---------------------------
		
		 
		% DMC, delU(k) �����ú�y1(k+1/k)
        DMC_y0 		= DMC_y0 	 + duk * DMC_step; 			% ����� y1
		y0p1		= y0p1       + duk * y_dmc_step(m_nHP+1);
		
		DMC_yideal	= DMC_yideal + duk * y_dmc_step;		% ����� y1,HP�㹻��
		
		%----  ������Ϣ	---------------------------
		DMC_y01_mem(itdmc,:) = DMC_y0';			% �洢DMC�������
		%----  ������Ϣ	---------------------------
		
    end
        % ���浽 �����У����ں���Ļ�ͼ
        d_dv(it)	= dk;		y_dv(it)	= ydk;		% �Ŷ�ͨ�� gd ���������
        u_mv(it)	= uk;    	y_pv(it)	= yk;		% ���ض��� gp ���������
        y_sum(it)	= yksum;

        t_time(it)	= tk;
end

figure();
msg = sprintf('û��ǰ��������DMC Controler: HP=%d, HM=%d, A(%d,1)=%f.LRJ����42', m_nHP, m_nHM, m_nHP,A(m_nHP,1));

% t_time = 0: h_gp : (itEnd-1)*h_gp;
	
subplot(211);
%plot(t_time, y_pv);
plot(t_time, y_pv, 'r-', t_time, y_dv, 'g-', t_time, y_sum, 'k-' );
ylabel('y');
legend('���ض�������ypv', '�Ŷ����yd', 'DMC�ı�����(ypv+yd)');

title(msg);
grid on;

subplot(212);
stairs(t_time, u_mv, 'r-');	hold on;
% stairs(t_time, d_dv,'g-');
legend('DMC��������u_mv', '�Ŷ����� d_dv');
ylabel('u');
xlabel('time');
title('�������� u, �Ŷ����� d, dead=0');
grid on;