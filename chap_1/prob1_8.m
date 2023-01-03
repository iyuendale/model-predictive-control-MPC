%% 1
num = 1; den = [1 2*0.001*1 1];
systf = tf(num, den)
[Am, Bm, Cm, Dm] = tf2ss(num, den)
% Am = [0 1; -1 -0.002], Bm= [0 0.1]', Cm = [0 1], Dm = 0
[Am, Bm, Cm, Dm] = c2dm(Am, Bm, Cm, Dm, .5), Bd = [0 1]';
A = [Am zeros(size(Am, 1), 1); Cm*Am 1]
B = [Bm; Cm*Bm], C = [zeros(1, size(Cm, 2)) 1]
%% 2
Np = 60; Nc = 20; rw = 0.1; R = rw*eye(Nc, Nc);
rs = 1; Rs = ones(Np, 1);
% observer design
poles = [0.1 0.2 0.3];
Kob = place(A', C', poles)'
%% MPC parameters
F= []; phi= [];
for k = 1:Np
	phi_col = [];
	for j= 1:Nc
		if k<j
			phi_col= [phi_col zeros(size(C, 1), size(B, 2))];
		else
			phi_col= [phi_col C*A^(k-j)*B];
		end
	end
	phi= [phi; phi_col];
	F = [F; C*A^k];
end

Ky = inv(phi'*phi+R)*phi'*Rs; Kmpc = inv(phi'*phi+R)*phi'*F;
x_old= [0 0]'; xe_old = [0 0.1 0.10]'; u = 0; d= 0;
xe = [0.1 0.1 0.1]'; x = [0 0 0]'; xm = [0 0]'
buf = [0 xe(3) x(3) 0 xm(2)]; bufu = [];
deltax = [0 0 0]';
for k = 1:50
	deltaU = Ky - Kmpc*xe;
	xe = A*xe + B*deltaU(1) + Kob*C*(x-xe);
	x = A*x + B*deltaU(1);
	u = u+deltaU(1);
	xm = Am*xm + Bm*u;
	buf = [buf; k xe(3) x(3) deltaU(1) xm(2)];
	bufu = [bufu; [k k+1]' [u u]'];
end
% k = buf(:, 1); X = buf(:, 3); Xe = buf(:, 2); deltaU = buf(:, 4);
% plot(k, [X, Xe]), hold on
% % for k = 100:120
% 	deltaU = Ky - Kmpc*xe;
% 	xe = A*xe + B*deltaU(1) + Kob*C*(x-xe);
% 	x = A*x + B*deltaU(1) + [0 0 1]'
% 	buf = [buf; k xe(3) x(3) deltaU(1)];
% end
k = buf(:, 1); x = buf(:, 3); xe = buf(:, 2); deltaU = buf(:, 4);
xm = buf(:, 5); u = bufu(:, 2); ku = bufu(:, 1);
subplot 411, plot(k, x)
subplot 412, plot(k, xe)
subplot 413, plot(k, xm)
subplot 414, plot(ku, u)


% 
% deltaU = Ky - Kmpc*xe_old;
% 	u = u+deltaU(1);
% 	xe_new = A*xe_old + B*deltaU(1) + Kob*C*(deltax-xe_old);
% % 	deltaxe = [(xe_new-xe_old)' ; xe_new(2)];
% 	x_new = Am*x_old + Bm*u + Bd*d;
% 	deltax = [(x_new-x_old); x_new(2)]
% 	buf = [buf; k xe_new(3) x_new(2) deltaU(1) u];
% 	x_old = x_new; xe_old = xe_new;