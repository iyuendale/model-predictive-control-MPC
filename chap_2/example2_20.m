%% discretize
am= [1.9048 -0.9048 0.0048; 1 0 0; 0 0 0];
bm= [0.0048 0 1]';
cm= [1 0 0];
dm= 0;

%% augmented system
A= [am [0 0 0]'; cm*am 1]; B= [bm; cm*bm];
C= [zeros(size(cm, 1), size(am, 1)), 1]; D= 0;
% eig(A)
%% F and phi
Np = 46; Nc= 8; R= 0.1*eye(Nc);
[F, phi]= F_phi(A, B, C, Np, Nc);
Rs= 1*ones(Np, 1);
phi'*phi; phi'*F;
E= phi'*phi + R;

%% with constraints applied on first element of Y
% M = [-C*B zeros(1, 7); C*B zeros(1, 7)];
M = [[-1 zeros(1, Np-1)]*phi; [1 zeros(1, Np-1)]*phi];
x= [0 0 0 0]'; buf = []; u = 0; buf2= []; xm = [0 0 0]'; xm_old= xm;
for k = 1:19
	FF= phi'*(F*x - Rs);
	b = [[1 zeros(1, Np-1)]*F*x; 1 - [1 zeros(1, Np-1)]*F*x];
	deltaU = QPhild(E, FF, M, b);
	u = u+deltaU(1);
	xm = am*xm + bm*u;
	x = [xm-xm_old; cm*xm];
	xm_old= xm;
	buf = [buf; k, cm*xm deltaU(1) u];
	buf2= [buf2; [k k+1]' [u u]'];
end
% input disturbance
for k = 20:39
	FF= phi'*(F*x - Rs);
	b = [1 zeros(1, Np-1); -1 zeros(1, Np-1)]*F*x + [0 1]';
	deltaU = QPhild(E, FF, M, b);
	u = u+deltaU(1);
	xm = am*xm + bm*u + [0 0 1]';    % a positive disturbance
	x = [xm-xm_old; cm*xm];
	xm_old= xm;
	buf = [buf; k, cm*xm deltaU(1) u];
	buf2= [buf2; [k k+1]' [u u]'];
end
% 
for k = 40:60
	FF= phi'*(F*x - Rs);
	b = [1 zeros(1, Np-1); -1 zeros(1, Np-1)]*F*x + [0 1]';
	deltaU = QPhild(E, FF, M, b);
	u = u+deltaU(1);
	xm = am*xm + bm*u - [0 0 1]';      % a negative disturbance
	x = [xm-xm_old; cm*xm];
	xm_old= xm;
	buf = [buf; k, cm*xm deltaU(1) u];
	buf2= [buf2; [k k+1]' [u u]'];
end
% 
k = buf(:, 1); y = buf(:, 2); deltaU = buf(:, 3); U = buf(:, 4);
k2 = buf2(:, 1); u = buf2(:, 2); xm = [0 0]';
subplot 311, plot(k, y), hold on % , axis([0 60, 0 1.5]), hold on
subplot 312, plot(k, deltaU), hold on % , axis([0 60 -2 7]), hold on
subplot 313, plot(k, U), hold on % , axis([0 60 -6 15]), hold on
% %% without constraints
% xm= [0 0 0]'; xm_old= xm; x= [0 0 0 0]'; buf = []; u = 0; buf2= [];
% 
% for k = 1:19
% 	FF= phi'*(F*x - Rs);
% 	deltaU= -E\FF;
% 	u = u+deltaU(1);
% 	xm = am*xm + bm*u;
% 	x = [xm - xm_old; cm*xm];
% 	buf = [buf; k, cm*xm deltaU(1) u];
% 	buf2= [buf2; [k k+1]' [u u]'];
% end
% 
% for k = 20:39
% 	FF= phi'*(F*x - Rs);
% 	deltaU= -E\FF;
% 	u = u+deltaU(1);
% 	xm = am*xm + bm*u + [0 0 1]';
% 	x = [xm - xm_old; cm*xm];
% 	buf = [buf; k, cm*xm deltaU(1) u];
% 	buf2= [buf2; [k k+1]' [u u]'];
% end
% 
% for k = 40:60
% 	FF= phi'*(F*x - Rs);
% 	deltaU= -E\FF;
% 	u = u+deltaU(1);
% 	xm = am*xm + bm*u - [0 0 1]';
% 	x = [xm - xm_old; cm*xm];
% 	buf = [buf; k, cm*xm deltaU(1) u];
% 	buf2= [buf2; [k k+1]' [u u]'];
% end

k = buf(:, 1); y = buf(:, 2); deltaU = buf(:, 3);
k2 = buf2(:, 1); u = buf2(:, 2);
subplot 311, plot(k, y), axis([0 60, 0 1.5]), title Output
legend 'constraints on first elements' 'constraints on all elements'
subplot 312, plot(k, deltaU), title \DeltaU
legend 'constraints on first elements' 'constraints on all elements'
subplot 313, plot(k2, u), title 'control signal u'
legend 'constraints on first elements' 'constraints on all elements'