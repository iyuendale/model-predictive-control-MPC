num= 10; den= [1 0.1 3];
[ac, bc, cc, dc]= tf2ss(num, den);
%% discretize
[am, bm, cm, dm]= c2dm(ac, bc, cc, dc, 0.01);

%% augmented system
A= [am [0 0]'; cm*am 1]; B= [bm; cm*bm];
C= [zeros(size(cm, 1), size(am, 1)), 1]; D= 0;
% eig(A)
%% F and phi
Np = 20; Nc= 3; R= 0.01*eye(Nc);
[F, phi]= F_phi(A, B, C, Np, Nc);
Rs= 1*ones(Np, 1);
phi'*phi; phi'*F;
E= phi'*phi + R;

%% observer design
poles= [0.001 0.0015 0.002];
Kob= place(A' , C', poles)';

%% with constraints applied on all of elements of U
M = [1 0 0; 1 1 0; 1 1 1; -1 0 0; -1 -1 0; -1 -1 -1];
xest= [0 0 0]'; x= [0 0 0]'; buf = []; u = 0; buf2= []; xm = [0 0]';
buf3 = [];
for k = 1:60
	FF= phi'*(F*xest - Rs);
	b = [6-u; 6-u; 6-u; 3+u; 3+u; 3+u];
	deltaU = QPhild(E, FF, M, b);
	u = u+deltaU(1);
	xest= (A - Kob*C)*xest + B*deltaU(1) + Kob*cm*xm;
% 	x = A*x + B*deltaU(1);
	xm = am*xm + bm*u;
	buf = [buf; k, cm*xm deltaU(1)];
	buf2= [buf2; [k k+1]' [u u]'];
	buf3= [buf3; k, deltaU'];
end
k = buf(:, 1); y = buf(:, 2); deltaU = buf(:, 3);
k2 = buf2(:, 1); u = buf2(:, 2); xm = [0 0]';
figure(1)
subplot 311, plot(k, y), axis([0 60, 0 1.5]), hold on
subplot 312, plot(k, deltaU), axis([0 60 -2 7]), hold on
subplot 313, plot(k2, u), axis([0 60 -6 15]), hold on

figure(2), plot(buf3(:, 1), buf3(:, 2:4))
legend 'Deltau(k)' 'Deltau(k+1)' 'Deltau(k+2)'
title 'control changes'
%% with constraints on only te first element of U
M = [1 0 0; -1 0 0];
xest= [0 0 0]'; x= [0 0 0]'; buf = []; u = 0; buf2= [];
for k = 1:60
	b = [6-u; 3+u];
	FF= phi'*(F*xest - Rs);
	deltaU = QPhild(E, FF, M, b);
	u = u+deltaU(1);
	xest= (A - Kob*C)*xest + B*deltaU(1) + Kob*cm*xm;
% 	x = A*x + B*deltaU(1);
	xm = am*xm + bm*u;
	buf = [buf; k, cm*xm deltaU(1)];
	buf2= [buf2; [k k+1]' [u u]'];
end
k = buf(:, 1); y = buf(:, 2); deltaU = buf(:, 3);
k2 = buf2(:, 1); u = buf2(:, 2);
figure(1)
subplot 311, plot(k, y), axis([0 60, 0 1.5]), title Output
legend 'constraints on all elements' 'constraint on first element'
subplot 312, plot(k, deltaU), title 'DeltaU'
legend 'constraints on all elements' 'constraint on first element'
subplot 313, plot(k2, u), title 'control signal u'
legend 'constraints on all elements' 'constraint on first element'