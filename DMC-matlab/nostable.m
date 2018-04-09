
% ��һ���֣����ض���ķ���
% ���ض������
T0_gp=25;	k0_gp = 10;

h_gp = T0_gp / 50;			% ���ض��� ���� ʱ�䲽��
a_gp = exp(-h_gp / T0_gp); b_gp = k0_gp*(1-a_gp);
Tend =  10 * T0_gp / h_gp;

y_pv = zeros(Tend, 1);			% �洢������,��hΪ���沽��
ysp = 1;

uk = 0;
y_pv(1) = 0;                % ��ʼ״̬ y(0)=0
yss = 0.99* ysp;



for ( it = 1 : Tend )
	y_pv(it+1) = a_gp * y_pv(it) + b_gp * uk;
	
	ek = ysp - y_pv(it+1);
	uk = 10 * ek;
		
    end

    t_time = 0: h_gp: (Tend)*h_gp;
    figure();
    % plot(t_time, y_pv(1:itEnd), 'b-'); grid on; hold on;
	plot( y_pv ) ;
