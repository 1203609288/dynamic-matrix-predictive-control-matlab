
Ts = 10;
Tr = 2;		Tr2 = 100;

% h  = Ts / Tr;

a = exp(-Ts/Tr);	
a = exp(-Ts/Tr2);
a1 = a;				

u = 2 ;
y0 = 0.5;

y(1) 		= y0;
y_mphc(1) 	= y0;

for( it = 1 : 100 )

	y(it+1) = a * y(it) + (1-a) * u; 	% 仿真算法
	
end

	
	
for( it = 1 : 100 )
	y_mphc(it+1) = a1 * y0 + (1-a1) * u; 	% MPHC 的算法
	a1 = a1*a;	
end

t = 1: 101;

figure();
plot( t, y, 'r', t, y_mphc, 'b');
legend('' , '' );
grid;
