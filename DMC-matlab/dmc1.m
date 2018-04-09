clear;
clc;
% ��һ���֣����ض���ķ���
% ���ض������
T0_gp=25;	k0_gp = 10;

h_gp = T0_gp / 50;			% ���ض��� ���� ʱ�䲽�� 0.5
a_gp = exp(-h_gp / T0_gp); b_gp = k0_gp*(1-a_gp);

Ts = 2;
m_nHP = 10;
m_nHM = 3;

[ADMC, A, y_dmc_step, m_nHP, m_nHM ] = CalcStepResponseADMC( T0_gp, k0_gp, Ts, m_nHP, m_nHM , false);
m_nHP0 = length(y_dmc_step) ;

%load dmc1_log ADMC y_dmc_step
%load dmc1_d0_log
%dModel=d_dmc_step(1:m_nHP);

% �ڶ����֣��Ŷ�ͨ��
T0_gd=10;	k0_gd = 5;

h_gd = T0_gd / 40;			% �Ŷ�ͨ�� ���� ʱ�䲽���� 0.25
a_gd = exp(-h_gd / T0_gd); b_gd = k0_gd*(1-a_gd);



% �������֣����ض����DMC����, 

% ���ض�����������
ysp = 1;	% �趨�����ź�
uk = 0 ;    % mv
yk = 0      % ���ض������������

dk=0;		% �Ŷ��źţ����Ŷ����ڵ������
ydk=0;		% �Ŷ����ڵ����

itEnd = 500;
t_time =  zeros(itEnd,1);

u_mv =  zeros(itEnd,1);			% �洢��ʷ���ݣ�������ͼ
y_pv =  zeros(itEnd,1);			
y_d  = zeros(itEnd,1);
y_sum = zeros(itEnd,1);

DMC_y1 = zeros( size(y_dmc_step) );	% �����Ծ��̬����Ӧ,HP�㹻��

DMC_y0 = zeros( m_nHP, 1 );		% �� deltalU(k) ����ʱ��������Ӧ
DMC_ye = zeros( m_nHP, 1 );		% ƫ��
DMC_du = zeros( m_nHM, 1 );		% mv����
DMC_yw = ysp * ones( m_nHP, 1 );		% �ο��켣. ������򵥴���=1

DMC_du_mem = zeros(itEnd, m_nHM);	%DMC������̵Ĵ洢
DMC_y0_mem = zeros(itEnd, m_nHP);
DMC_y1_mem = zeros(itEnd, m_nHP0);

it = 1; itdmc=0;
yerror = 0;

clc;
while (it < itEnd)
	it = it + 1;
	
	%----	 �Ŷ����ڵļ���		----
	
	%	 �Ŷ������ź�
	if( 50<it && it<70)     dk = 5;		% �����ʩ���Ŷ������� dk = 0 
    else                    dk = 0;
	end
	
	% �Ŷ��Ĳ��� h_gd = 0.25, �� ���ض���ķ��沽�� h_gp = 0.5 
	% �����Ҫ����2��
	ydk = a_gd * ydk + b_gp * dk;
	ydk = a_gd * ydk + b_gp * dk;
    
	ydk = 0;
	
	%---- ���ض���ķ������	----
	yk = a_gp * yk + b_gp * uk ;

    if( yk <= 0.9 & it > 10 )
        it, uk, yk;
    end
	
	% y_sum = ���ض������� y_pv + �Ŷ����ڵ���� ydk
	yksum = yk + ydk;
	
	if( mod(it, 4) == 1 )	
        
        itdmc = itdmc + 1;
        
		if( itdmc == 1 )	% ��һ��ִ��DMC����
			DMC_y0 = yk*ones(size(DMC_y0));
			DMC_y1 = yk*ones(size(DMC_y1));
		else
			yerror = yksum - DMC_y0(1) ;			% ƫ��
			DMC_y0 ( 2 : m_nHP )   = DMC_y0 ( 2 : m_nHP ) + yerror; 	% ����У��				
			DMC_y0 ( 1 : m_nHP-1 ) = DMC_y0 ( 2 : m_nHP );			% ����
				
			% DMC_y0 (m_nHP) = DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2)+DMC_y0(m_nHP); % LRJ�Ľ���

			DMC_y1 ( 2 : m_nHP0 )   = DMC_y1 ( 2 : m_nHP0 ) + yerror; 
			DMC_y1 ( 1 : m_nHP0-1 ) = DMC_y1 ( 2 : m_nHP0 );
		end
		
		% DMC_yd = dModel * dk;				% ���Ŷ�ydk��һ��Ԥ��
		
		% �趨ֵ-Ԥ��y-ǰ������
		% DMC_ye = DMC_yw - DMC_y0 - DMC_yd;	% �趨ֵ-Ԥ��y-ǰ������

		% ��ʽû��ǰ������
		DMC_ye = DMC_yw - DMC_y0 ;		% �趨ֵ-Ԥ��y		
		
		DMC_du = ADMC * DMC_ye ;		% ������ƹ��� ��uk			
		uk = uk + DMC_du(1);
		
		% PID ����,  ��K>2ʱ, �������ֵ���㲻�ȶ�������,���·�������,��ɢ.
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
		 
        DMC_y0 = DMC_y0 + DMC_du(1) * A(1:m_nHP,1); 	% ����� y1
		DMC_y1 = DMC_y1 + DMC_du(1) * y_dmc_step;		% ����� y1,HP�㹻��      
		
    end
        % ���浽 �����У����� ����Ļ�ͼ
        u_mv(it) = uk;    
        y_d (it) = ydk;
        y_pv(it) = yk;	
        y_sum(it)= yksum;    
end

figure();
msg = sprintf('û��ǰ��������DMC Controler: HP=%d, HM=%d, A(%d,1)=%f', m_nHP, m_nHM, m_nHP,A(m_nHP,1));

subplot(211);
t_time = 0: h_gp : (itEnd-1)*h_gp;
%plot(t_time, y_pv);
plot(t_time, y_pv, 'r-', t_time, y_d, 'g-', t_time, y_sum, 'k-' );
ylabel('y');
title(msg);
legend('���ض�������ypv', '�Ŷ����ykd', 'DMC�ı�����(ypv+ykd)');
grid on;

subplot(212);
stairs(t_time, u_mv);
ylabel('u');
xlabel('time');
title('�������� uk');
grid on;
