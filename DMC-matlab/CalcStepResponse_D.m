% �Ŷ�ͨ�� ���棬��ȡ ��Ծ��Ӧ ����
% y_step     = ���ض��� �� ��̬ �������Ľ�Ծ��Ӧ
% y_dmc_step = y_dmc_stepy( 1 : m_nHP )
function [ y_dmc_step, y_step ] =  CalcStepResponse_D( T0_gd, k0_gd, Ts, m_nHP, m_nHM , bFigure)


disp('�Ŷ�ͨ���ķ���');

% ��һ���֣��Ŷ�ͨ���ķ���

% �Ŷ�ͨ������
% T0_gd=10;		k0_gd = 5;
h_gd = T0_gd / 50;			% ���� ʱ�䲽��,  0.2
a_gd = exp(-h_gd / T0_gd); b_gd = k0_gd*(1-a_gd);
iTend = 10 * T0_gd / h_gd;

y_pv		= zeros(iTend, 1);		% �洢������,��hΪ���沽��
t_time		= zeros(iTend, 1);

y_step		= zeros(m_nHP, 1);	% ��������,�� Ts Ϊ��������
y_dmc_step	= zeros(m_nHP, 1);	% ��������,�� Ts Ϊ��������
uk = 1;

for ( it = 1 : iTend-1 )	% while ( it < Tend )
	y_pv(it+1) = a_gd * y_pv(it) + b_gd * uk;	% y_pv(1) = t0, y_pv(2) = t1

	% t_time(it+1) = t_time(it) + h_gd;			% ��ֵ�����������һЩ
	t_time(it+1) = it * h_gd;

	%if( y_pv(it+1) > 0.98*k0_gd )
     %   itEnd = it+1;
     %  break;
    % end
end

itEnd = it+1;
t_time2 = [ 0: h_gd: (itEnd-1)*h_gd ]';				% ֱ������ʱ�����꣬�����ٶȸ���
% assert( t_time == t_time2 );						% ��ֵ����������²���ȣ�
if( bFigure)		% �����Ҫ����ͼ���� �Ŷ��Ľ�Ծ��Ӧ �������ģ�
	figure();
	plot(t_time, y_pv(1:itEnd), 'r-'); 		hold on;
	end

% �ڶ����� ���ݿ������Ĳ������ڣ����� �������õ� ��Ծ��Ӧ
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

y_dmc_step = y_step(1 : m_nHP);			% ֻҪ��Ŀ���ʱ�������
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

	msg = sprintf('�Ŷ�����=%g/(%gs+1), h=%g, m_nHP=%d', k0_gd, T0_gd, h_gd, m_nHP );
	title(msg,'Interpreter', 'none');


end

save CalcStepResponse_D y_dmc_step y_step