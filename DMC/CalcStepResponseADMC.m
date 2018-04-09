function [ADMC, A, y_step, m_nHP, m_nHM ] =  CalcStepResponseADMC( T0_gp, k0_gp, Ts, m_nHP, m_nHM , bFigure, nLoop, Kp2)

% ���ض��� ���棬��ȡ ��Ծ��Ӧ ����
% y_step = ���ض��� �� ��̬ �������Ľ�Ծ��Ӧ
% A		: ��̬����
% ADMC	: deltaU = ADMC * DMC_ye
% m_nHP	: Ԥ��ʱ�򳤶�
% m_nHM	: ����ʱ�򳤶�

%m_nHP = 100;
%m_nHM = 5;
%bFigure = true;

% ��һ���֣����ض���ķ���
% ���ض������
% 
% ��������ȱʡ����
% T0_gp=25;	k0_gp = 10;
% h_gp = T0_gp / 50;			% ���ض��� ���� ʱ�䲽��
% Tend =  10 * T0_gp / h_gp;

%---- 2013 �������¶���
%nLoop = 1;		% �ڻ�·
%nLoop = 2;		% ���·

n_gp=5;							% 5�׶���, 2�����ڻ�·��3�������·
%T0_gp = [15, 15,	20,     20, 20];
%k0_gp = [8,  1,		1.125,  1,  1 ];
h_gp  = 0.5;					% ���ض��� ���� ʱ�䲽��

a_gp = exp( - h_gp ./ T0_gp );
b_gp = k0_gp .* ( 1 - a_gp );
iTimeSimuEnd = 4500;			% ���泤��,
iTimeSS = 0;					% �ﵽ��̬, static state

% y_pv(:, [�ڻ�·��ǰ�����£� ���·����������] )
y_pv = zeros(iTimeSimuEnd, 2);	% �洢������,��hΪ���沽��
y_step = zeros(m_nHP,1);	    % ��������,�� Ts Ϊ��������
y_pv(1,1) = 0;                % ��ʼ״̬ y1(0)=0
y_pv(1,2) = 0;                % ��ʼ״̬ y2(0)=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%�����������ֵ%%%%%%%%%%%%%%%%%%%%%%%%%%%

uk = 1;
if ( nLoop ==1 )		% ֻ�ڻ�· ��
	yss = 0.9999 * k0_gp(1) *k0_gp(2) * uk;		% ��ǰ������ֵ̬
	nypvCol = 1;
elseif( nLoop == 2 )	% ֻ���· ��
	yss = 0.9999 * k0_gp(1) *k0_gp(2) * k0_gp(3) *k0_gp(4) * k0_gp(5) * uk;	% ���ڶ���������
	nypvCol = 2;
elseif( nLoop == 3 )	% �ڻ�· P�� ���· DMC
	Kinss = ( k0_gp(1) * k0_gp(2) * Kp2 );	
	Kinss = Kinss / (1 + Kinss );
	% yss = 0.9999 * Kinss / (1 + Kinss ) * uk ;	 % �� �ڻ�·����ֵ̬����Ϊ ss
	yss = 0.9999 * Kinss *  k0_gp(3) *k0_gp(4) * k0_gp(5) * uk;	
	nypvCol = 2;
end	

x = zeros(n_gp+1,1);		% �����м���� x
 x(1) = uk;
for ( it = 1 : iTimeSimuEnd-1 )	
	if( nLoop ==1 | nLoop== 2)			
		x(1) = uk;			% �ڻ�·����������
	elseif( nLoop == 3 )	
		x(1) = Kp2 *( uk - x(3) ) ;		% �ڻ�·����P���ƣ�KP=Kp2
	end
	
	for( inh = 1 : 5 )
		x(inh+1) = a_gp(inh) * x(inh+1) + b_gp(inh) * x(inh);
	end	
	y_pv(it+1,1) = x(3);	y_pv(it+1,2) = x(6);
	
%	�ﵽ��̬�����˳�����
%	if( y_pv(it+1,1) > yss2 )
%       iTimeSS = it+1;
%        break;
%	end	
	if( (y_pv(it+1,nypvCol) > yss)  & (iTimeSS == 0 ) )
        iTimeSS = it+1;
        break;
	end	
end
y_len = max ( it+1, iTimeSS) ;	% ��ͼʱ������ȡ���������

if( bFigure )
    t_time = 0: h_gp: (y_len-1)*h_gp;
    figure();
    plot(t_time, y_pv(1:y_len,1), 'r-', t_time, y_pv(1:y_len, 2), 'k-');
	% legend('��ǰ�����£��ڻ�·', '���������£����·', 'Location', 'SouthEast');	
	grid on;
    hold on;
end


% �ڶ����֣����ݿ������Ĳ������ڣ����� �������õ� ��Ծ��Ӧ
% Ts = 2;			% �������� =2 
mulTime = Ts / h_gp;

it2 = 1;			% ��Ծ��Ӧ�� ѭ�� a(it2)
it  = 1+mulTime;	% 4 = Ts / h 
while( it <= iTimeSS )
    y_step(it2) = y_pv(it,nypvCol);		
	
    it = it + mulTime;
    it2 = it2 + 1;
end
% y_step ( it2 ) = y_pv(iTimeSS);     /// �Ƿ� aN = ass

m_nHP = min(m_nHP, it2-1);

% �������֣�����DMC���õ� ģ�;��� A 
A=zeros(m_nHP, m_nHM);	% DMC��ģ�;��� A 
for( it = 1 : m_nHM )
	A([it : m_nHP], it) = y_step( [1 : m_nHP+1-it] );
end

if( bFigure)
    t2_time = Ts : Ts : (it2-1)*Ts;
	plot(t2_time, y_step(1 : length(t2_time) ), 'b.');
    plot(t2_time([1:m_nHP]), A([1:m_nHP],1), 'm*');
    grid on;
	legend('��ǰ�����£��ڻ�·', '���������£����·', 'y.step(Ts)', 'y.dmc.step(m.nHP*Ts)', 'Location', 'SouthEast');
	
	xlabel('time');
	ylabel('y');
% 	msg = sprintf('���ض���: nLoop==%d, h=%g, Ts=%g, m_nHP=%d, m_nHM=%d',nLoop, h_gp, Ts, m_nHP, m_nHM );
% 	title(msg,'Interpreter', 'none');
title('���Ȳ��ֽ�Ծ��Ӧ�仯');
end

%if( m_nHP < length(y_step ) )
%    A(m_nHP, 1) = 10; % //y_step( length(y_step) );
%end

% ���Ĳ��֣�����DMC���õ� ģ�;��� A_DMC 

% DMC �Ŀ������Ŀ��Ʋ���
m_fQ = eye( m_nHP, m_nHP );
% for( i=1:m_nHP )     m_fQ(i,i)= i*i; end
m_fR = 5*eye( m_nHM, m_nHM );
% m_fR = zeros( m_nHM, m_nHM ) ;

% ���棬���� DMC ���� DMC = inv(A'QA+R) A'Q
At = A';				% A ��ת��'
At1 = At * m_fQ;
At2 = At1 * A;
At3 = At2 + m_fR;
At4 = inv(At3);
At5 = At4*At1;
ADMC = At5

%A(m_nHP, 1) = 10;
% save dmc1_log ADMC y_step

% y_step ( m_nHP+1 , 1 ) = uk * k0_gp ;                % = 10, �о�ass������

if( nLoop==1)
	save CalcStepResponseADMC_in  y_step A ADMC
else
	save CalcStepResponseADMC_out y_step A ADMC
end
	