% ���ļ���ʾ���������DMC�㷨	2011.08

% �� (k+1)ʱ�̣�����du(k+1)ʱ,  ��Ҫ y0(k+1+i/k+1), i=1,2,...HP.
%   i=1,2,...(HP-1),  y0(k+1+i/k+1) = y1(k+1+i/k) + h*e = y0(k+1+i/k) + a(1+i) * du(k) + h*e
%   i=HP,  ��ΪԤ��ʱ��=HP, ֻ�� y0(k+HP),  ��û�� y0(k+1+HP)  
%  	�� ����� y(HP+1)=y(HP)������, 
%   ��׼DMC�㷨�У�HPȡ�㹻����ʹ��ϵͳ�ﵽ��̬. �����Ϊ  y(HP+1)=y(HP)�� �� y0(k+1+HP/k+1) = y0(k+HP/k+1)
%   ���HP����������y0(k+1+HP/k+1) ������ y0(k+HP/k+1)�� �ڵ��ڹ��̽���ʱ���ᷢ��΢С���񵴡�Ϊ����������,	LRJ��ʾ����������������

% LRJ��������1:	��Ϊ ��Ծ��Ӧ��������ƽ����,  ����ʽ�� ���ۣ�����������Ϊ�ǳɹ���,	ʹ�� ϵͳƽ�� HP �����ں�Ĳ�����С
%  DMC_y0 (m_nHP) += DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2) %

% LRJ��������2:	 �������������ۣ������Ͽ�ȡ���ɹ�. �����������̫��, ���Ʊ�ʹ��ϵͳ�Ŀ��������ܽ���
		%if( -0.6<duk &	duk<0.6	)
		%	duk	= 0	;
		%end

% LRJ��������3:	���ɹ�
%	DMC_step(m_nHP)	= y_step(m_nHP-1);		% ������ y1(k+HP/k)	ʱ�õĽ�Ծ��Ӧϵ�� �޸���.
%	��Ϊ�ڼ��� y1(k+1+HP/k+1) =	y0(k+1+HP/k+1) + a(HP) * du(k+1) �������
%	�ʲ���	   y1(k+1+HP/k+1) =	y0(k+1+HP/k+1) + a(HP-1) * du(k+1),	�ǿ��ǵ� ��ʵ  y0(k+1+HP/k+1)����ֵ�� =	 y0(k+HP/k+1)
%	������Ϊ������ԭ��,	y0(k+1+HP/k+1) = y0(k+HP/k+1), ���	y1(k+1+HP/k+1) ʵ������ֵ��= y1(k+HP/k+1)
%
% LRJ��������4:	���ǵ�������������,Ҳ�ǲ��ɹ���.
%	y0(HP+1) �ĳ���Ϊ HP+1,
%	��ʹ��DMC���Ƶ�ʱ��, ʱ��Ϊ	HP.
%	������������ֻ�ǰ�С�����Ӻ���һ������,����ʵ���Եĸı�. �൱�� HP������,	����Ч����һЩ.
%
% ���ļ�,��ʾ��ǰ���ַ���
%
% ���ļ�,��ǰ������LRJ��������1

clear;
clc;
bShowMsg = false;
% ��һ���֣����ض���Ĳ���,���ڽ������ض����ģ��,�Լ�DMC���Ƶķ���.

%T0_gp=25;	k0_gp =	10;
%h_gp = T0_gp / 50;			% ���ض���������.	 ����ʱ�䲽��=0.5
%a_gp = exp(-h_gp / T0_gp); b_gp	= k0_gp*(1-a_gp);

% ���¶��󱻿ض������Gp
%
%			8				%			1.125
%	------------------		%	------------------
%		(1+15s)^2			%		(1+20s)^3

n_gp=75;							% 5�׶���, 2�����ڻ�·��3�������·
% T0_gp = [15 , 15 , 20 , 20 , 20, 0];
% k0_gp = [8,  1, 1.125 ,  1,   1 ,1];
% T0_gp = [33.6,33.6,38.2, 38.2,	38.2, 38.2]; %%����A��
% k0_gp = [2.6,  1, 2.6 ,  1,   1 ,1];
% T0_gp = [88.65 , 88.65 , 35.46 , 35.46 , 0];%һ��
% k0_gp = [0.5625, 1, 0.8639 ,  1, 1];
% T0_gp = [22.5 , 22.5 , 24.6 ,  24.6,   24.6 ,24.6,24.6];%%���ȶ���B��
% k0_gp = [3.8,  1, 3.3 , 1, 1 , 1, 1];
% T0_gp = [88.65 , 88.65 ,153.51 , 153.51, 0];
% k0_gp = [9.975,  1, 1.722 , 1, 1 ];
T0_gp = [38.5 , 38.5 ,84.5 , 84.5, 0];
k0_gp = [7.75,  1, 1.25 , 1, 1 ];


h_gp  = 0.5;					% ���ض��� ���� ʱ�䲽��
a_gp = exp( -h_gp ./ T0_gp );
b_gp = k0_gp .* ( 1 - a_gp );

Ts = 2;						% DMC�ջ�����,��������
m_nHP =	400;				% DMC ����,	Ԥ��ʱ��
m_nHM =	50;					% DMC ����,	����ʱ��

% ���ض���ķ��棬 ��ȡ��Ծ��Ӧ
nLoop =2 ;					% 1 �ڻ�·��2=���·
bShowFigure=true;
[ADMC, A, y_step, m_nHP, m_nHM ] = CalcStepResponseADMC( T0_gp,	k0_gp, Ts, m_nHP, m_nHM	, bShowFigure, nLoop);
% load CalcStepResponseADMC_in
% load CalcStepResponseADMC_out

m_nHPideal = length(y_step)	;	% ����Gp�ﵽ��̬, ������Ԥ��ʱ��,

DMC_step = y_step(1:m_nHP);
% DMC_step(m_nHP) =	y_step(m_nHP-1);	% LRJ��������3:	���ɹ�


% �ڶ����֣��Ŷ�ͨ��, ���¶�����ʱ�������Ŷ�
b_FeedForward = false;
T0_gd=10;	k0_gd =	5;
h_gd = T0_gd / 40;				% �Ŷ�ͨ������.	����ʱ�䲽��=0.25
a_gd = exp(-h_gd / T0_gd); b_gd	= k0_gd*(1-a_gd);

% �Ŷ�ͨ���ķ��棬 ��ȡ��Ծ��Ӧ
%[DMC_D_step, d_step	] =	CalcStepResponse_D(T0_gd, k0_gd, Ts, m_nHP,	m_nHM, true);

% �������֣����ض����DMC����,

% ���ض�����������
ysp	= 1;	% �ջ�����,	�趨�����ź�
duk	= 0;	% u(k) - u(k-1)
uk = 0 ;	% mv
yk = 0 ;	% ���ض�������
y2k = 0;	% �ڻ�·(��ǰ��)���ض�������

dk=0;		% �Ŷ��źţ����Ŷ����ڵ������
ydk=0;		% �Ŷ����ڵ����

iTimeSimuEnd =  4500;			% ���泤��,
%itEnd	= 500;					% ���泤��
t_time	= zeros(iTimeSimuEnd,1);		% ����ʱ������,	��0��ʼ. t_time(1)=0, t_time(2)=h_gp, t_time(iTimeSimuEnd)	= (iTimeSimuEnd-1)*h_gp

u_mv	= zeros(iTimeSimuEnd,1);		% uk, �洢��ʷ����,������ͼ.
y_pv	= zeros(iTimeSimuEnd,1);		% yk, 
y2_pv 	= zeros(iTimeSimuEnd,1);		% yk, ��ǰ����, 

d_dv	= zeros(iTimeSimuEnd,1);		% dk
y_dv	= zeros(iTimeSimuEnd,1);		% ydk

y_sum	= zeros(iTimeSimuEnd,1);		% yk + ydk,

DMC_yideal	= zeros( size(y_step) );	% �����Ծ��̬����Ӧ, HP�㹻��,	����״̬. ������������ DMC_y0 �Ա�

DMC_y0	= zeros( m_nHP,	1 );			% �� deltaU(k) ����ʱ��������Ӧ, y0(k+i/k).		�����deltaU(k)��,�򱻸���,�൱�� y1(k+i/k)	= y0(k+i/k+1)
DMC_ye	= zeros( m_nHP,	1 );			% ƫ��,	w(k+i/k) - y(k+i/k)
DMC_du	= zeros( m_nHM,	1 );			% deltaU(k), mv����
DMC_yw	= ysp *	ones( m_nHP, 1 );		% �ο��켣.	w(k+i/k),	������򵥴���=1

DMC_du_mem = zeros(iTimeSimuEnd, m_nHM);		%DMC������̵Ĵ洢,���ڷ����������
DMC_y0_mem = zeros(iTimeSimuEnd, m_nHP);
DMC_y01_mem	= zeros(iTimeSimuEnd, m_nHP);		%�����deltaU(k)��,�򱻸���,�൱�� y1(k+i/k) = y0(k+i/k+1)
DMC_ye_mem = zeros(iTimeSimuEnd, m_nHP);
DMC_er_mem = zeros(iTimeSimuEnd, m_nHP);		% y(k) - y(k/k-1)

DMC_yideal_mem = zeros(iTimeSimuEnd, m_nHPideal);

it = 1;	itdmc=0;				% ѭ������
tk = 0;
yerror = 0;						% �ջ�У��,	ʵ����� = y(k+1) -	y1(k+1/k)

clc;

x = zeros(n_gp+1,1);		% �����м���� x
while (it <	iTimeSimuEnd)
	it = it	+ 1;				% ��һ��, it=2,	y(2)=t1ʱ��, y(1)=0ʱ�̵ĳ�ֵ.
	tk = tk	+ h_gp;				% ��ǰʱ��,	���ڵ���

	%----	 �Ŷ����ڵļ���		----		
	if(	b_FeedForward )
		%	 �Ŷ������ź�
		if(	100<it && it<200)	dk = 0.2;		% �����ʩ���Ŷ������� dk =	0
		else					dk = 0;
		end
		
		% �Ŷ��Ĳ��� h_gd =	0.25, �� ���ض���ķ��沽��	h_gp = 0.5,	�����Ҫ��������
		ydk	= a_gd * ydk + b_gp	* dk;
		ydk	= a_gd * ydk + b_gp	* dk;
	else
		ydk =	0;					% ����������Ŷ��ź�,���� =0
	end

	%---- ���ض���ķ������	----
	%yk = a_gp *	yk + b_gp *	uk ;
	x(1) = uk;	
	% x(1) = 2.5678 *(uk - x(3)) ;	
	% x(1) = 5 *(uk - x(3));	
	for( inh = 1 : 5 )
		x(inh+1) = a_gp(inh) * x(inh+1) + b_gp(inh) * x(inh);
	end
	
	y2k = x(3);		yk  = x(6);
	
	%if( yk	<= 0.9 & it	> 10 )
	 %	 it, uk, yk;
	%end
	
	yksum =	yk + ydk;			% yksum	= ���ض������� y_pv +	�Ŷ����ڵ���� ydk

	if(	mod(it,	4) == 1	)

		itdmc =	itdmc +	1;		% ��һ��, itdmc=1

		if(	itdmc == 1 )		% ��һ��ִ��DMC����
			DMC_y0 = yk* ones( size(DMC_y0)	);
			DMC_yideal = yk* ones( size(DMC_yideal)	);
		else
			yerror = yksum - DMC_y0(1) ;								% ʵ��ƫ��
			DMC_y0 ( 2 : m_nHP )   = DMC_y0	( 2	: m_nHP	) +	yerror;		% ����У��
			DMC_y0 ( 1 : m_nHP-1 ) = DMC_y0	( 2	: m_nHP	);				% ����
			DMC_er_mem(itdmc,1)	= yerror;
			% LRJ��������1:	��Ϊ ��Ծ��Ӧ����, ������Ϊ�ǳɹ���
			DMC_y0 (m_nHP) = DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2)+DMC_y0(m_nHP);	%
																																																																																																																																													  %	DMC_y0 (m_nHP) = DMC_y0(m_nHP-1)-DMC_y0(m_nHP-2) + DMC_y0(m_nHP);
			yerror		= yksum	- DMC_yideal(1);
			DMC_yideal ( 2 : m_nHPideal	  )	= DMC_yideal ( 2 : m_nHPideal )	+ yerror;
			DMC_yideal ( 1 : m_nHPideal-1 )	= DMC_yideal ( 2 : m_nHPideal );
		end
		
		if( b_FeedForward )							% ��ǰ������
			DMC_yd = DMC_D_step	* dk;					% ���Ŷ�ydk��һ��Ԥ��			
			% �趨ֵ - Ԥ��y  -ǰ������
			DMC_ye = DMC_yw -	DMC_y0 - DMC_yd;	% �趨ֵ-Ԥ��y-ǰ������
		else
			DMC_yd = 0;	
			% ��ʽû��ǰ������
			DMC_ye = DMC_yw	- DMC_y0 ;		% �趨ֵ-Ԥ��y
		end

		DMC_du = ADMC *	DMC_ye ;		% ������ƹ��� ��uk
		duk	= DMC_du(1);

		% LRJ��������2:	 ���������������Ͽ�ȡ���ɹ�
		%if( -0.6<duk &	duk<0.6	)
		%	duk	= 0	;
		%end
		uk = uk	+ duk;

		% PID ����,
		% ��K>2ʱ, �������ֵ���㲻�ȶ�������,���·�������,��ɢ.
		% uk = 10 *	(ysp - yksum);

		%----  ������Ϣ	---------------------------
		%----  ������Ϣ,����ƫ��
		if ( bShowMsg )
			msg	= sprintf('it=%d, itdmc= %d, time=%g, \t yerror=%f,	' ,	it,	itdmc, (it-1)*h_gp,	yerror)	;
			disp(msg);

			% [[1:50]',DMC_y0, DMC_ye];
			if(	DMC_du(1) >	0.01 &&	it > 10)
				[u_mv(it-10:it), y_pv(it-10:it)]
				DMC_ye,
			end
		end
		%----  ������Ϣ,����ƫ��

		% DMC������̵Ĵ洢,���ڷ����������
		DMC_ye_mem(itdmc,:)	= DMC_ye';
		DMC_y0_mem(itdmc,:)	= DMC_y0';
		DMC_yideal_mem(itdmc,:)	= DMC_yideal';

		DMC_du_mem(itdmc,:)	 = DMC_du';
		%----  ������Ϣ	---------------------------


		% DMC, delU(k) �����ú�y1(k+1/k)
		DMC_y0		= DMC_y0	 + duk * DMC_step;		% ����� y1

		DMC_yideal	= DMC_yideal + duk * y_step;		% ����� y1,HP�㹻��

		%----  ������Ϣ	---------------------------
		DMC_y01_mem(itdmc,:) = DMC_y0';			% �洢DMC�������'
		%----  ������Ϣ	---------------------------

	end
		% ���浽 �����У����ں���Ļ�ͼ
		d_dv(it)	= dk;		y_dv(it)	= ydk;		% �Ŷ�ͨ�� gd ���������
		u_mv(it)	= uk;		y_pv(it)	= yk;		% ���ض��� gp ���������
		y2_pv(it)	= y2k;		% �ڻ�·		
		y_sum(it)	= yksum;

		t_time(it)	= tk;
end

figure();
%msg = sprintf('û��ǰ��������DMC Controler: HP=%d,	HM=%d, A(%d,1)=%f',	m_nHP, m_nHM, m_nHP,A(m_nHP,1));
% msg	= sprintf('������DMC Controler:	HP=%d, HM=%d, A(%d,1)=%f��LRJ��������1. ����ΪCV.', m_nHP,	m_nHM, m_nHP,A(m_nHP,1));

subplot(211);
hl = plot(t_time, y2_pv, 'b-', t_time, y_pv, 'r-', t_time, y_dv, 'g-', t_time, y_sum,	'k-' );
set( hl(2), 'linewidth', 2.5);
legend('��ǰ��,�ڻ�·����ǰ����y2', '������,���·,������ypv, CV', '�Ŷ����yd', 'DMC�ı�����(ypv+yd)');

ylabel('y');
title('���Ȳ���DMC���Ʊ仯����');
grid on;


subplot(212);
stairs(t_time, u_mv, 'r-');	hold on;
stairs(t_time, d_dv,'g-');
legend('DMC��������u_mv', '�Ŷ�����	d_dv');
ylabel('u');
xlabel('time');
% title('��������	u, �Ŷ�����	d	deltaU��������0');
title('��������');
grid on;