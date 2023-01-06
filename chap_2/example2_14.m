num= 10; den= [1 0.1 3];
[ac, bc, cc, dc]= tf2ss(num, den);
%% discretize
[a, b, c, d]= c2dm(ac, bc, cc, dc, 0.01);

%% augmented system
A= [a [0 0]'; c*a 1]; B= [b; c*b];
C= [zeros(size(c, 1), size(a, 1)), 1]; D= 0;
eig(A)
%% F and phi
Np = 20; Nc= 3; R= 0.01*eye(Nc);
[F, phi]= F_phi(A, B, C, Np, Nc);
Rs= 1*ones(Np, 1);
phi'*phi; phi'*F;
E= phi'*phi + R;

%% observer design
poles= [0.001 0.0015 0.002];
Kob= place(A' , C', poles)';

%% with constraints applied on the first element of deltaU
M = [1 0 0; -1 0 0]; b = [3 1.5]';
xest= [0 0 0]'; x= [0 0 0]'; buf = []; u = 0; buf2= [];
for k = 1:60
	FF= phi'*(F*xest - Rs);
	deltaU = QPhild(E, FF, M, b);
	u = u+deltaU(1);
	xest= (A - Kob*C)*xest + B*deltaU(1) + Kob*C*x;
	x = A*x + B*deltaU(1);
	buf = [buf; k, x(3) deltaU(1)];
	buf2= [buf2; [k k+1]' [u u]'];
end
k = buf(:, 1); y = buf(:, 2); deltaU = buf(:, 3);
k2 = buf2(:, 1); u = buf2(:, 2);
subplot 311, plot(k, y), axis([0 60, 0 1.5]), hold on
subplot 312, plot(k, deltaU), axis([0 60 -2 7]), hold on
subplot 313, plot(k2, u), axis([0 60 -6 15]), hold on
%% without constraints
xest= [0 0 0]'; x= [0 0 0]'; buf = []; u = 0; buf2= [];
for k = 1:60
	FF= phi'*(F*xest - Rs);
	deltaU = -E\FF;
	u = u+deltaU(1);
	xest= (A - Kob*C)*xest + B*deltaU(1) + Kob*C*x;
	x = A*x + B*deltaU(1);
	buf = [buf; k, x(3) deltaU(1)];
	buf2= [buf2; [k k+1]' [u u]'];
end
k = buf(:, 1); y = buf(:, 2); deltaU = buf(:, 3);
k2 = buf2(:, 1); u = buf2(:, 2);
subplot 311, plot(k, y), axis([0 60, 0 1.5]), title Output
legend 'constrained' 'no constraint'
subplot 312, plot(k, deltaU), title \DeltaU
legend 'constrained' 'no constraint'
subplot 313, plot(k2, u), title 'control signal u'
legend 'constrained' 'no constraint'