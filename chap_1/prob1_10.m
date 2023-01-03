Am= [0.9048 0; 0.0952 1], Bm= [0.0952; 0.0048], Cm= [0 1], Dm = 0
%% 1
% Augmented system
A = [Am zeros(size(Am, 1), size(Cm,1)); Cm*Am eye(size(Cm, 1))]
B = [Bm; Cm*Bm], C = [zeros(size(Cm, 1), size(Am, 1)) eye(size(Cm, 1))]
D = 0
% controllability and observability checks
C_check = ctrb(A, B), C_rank = rank(C_check)
if C_rank== size(A, 1)
	disp('Augmented system is controllable')
else
	disp('Augmented system is uncontrollable')
end
O_check = obsv(A, C), O_rank= rank(O_check)
if O_rank== size(A)
	disp('Augmented system is observable')
else
	disp('Augmented system is unobservable')
end
%% 2
Np = 10; Nc = 4; R = 0.1*eye(Nc); 
F = []; phi = []; RRs = ones(Np, 1);
for k = 1:Np
	phi_row = [];
	for j = 1:Nc
		if k < j
			phi_row = [phi_row zeros(size(C, 1), size(B, 2))];
		else
			phi_row = [phi_row C*A^(k-j)*B];
		end
	end
	F= [F; C*A^k]; phi = [phi; phi_row];
end
F, phi
Ky = [1 zeros(1, Nc-1)]*inv(phi'*phi+R)*phi'*RRs
Kmpc = [1 zeros(1, Nc-1)]*inv(phi'*phi+R)*phi'*F
eig_closed= eig(A-B*Kmpc)
%% 3
poles = [0.1, 0.2, 0.3];
Kob= place(A', C', poles)'
%% 4
r = 2; d = 0;
x = [0 0]'; xe = [0 0 0]'; buf = []; u = 0;
bufu = [];
for k = 1:19
	deltaU = Ky*r - Kmpc*xe;
	xe = A*xe+B*deltaU+Kob*(x(2)-xe(3));
	u = u+deltaU;
	x = Am*x+Bm*u+[0 1]'*d;
	buf = [buf; k u deltaU xe(3) x(2)];
	bufu = [bufu; [k k+1]' [u u]'];
end
k = buf(:, 1); y = buf(:, 5); u = buf(:, 2);
% Disturbance rejection
d = -0.5; u = u(end);
for k = 20:60
	deltaU = Ky*r - Kmpc*xe;
	xe = A*xe+B*deltaU+Kob*(x(2)-xe(3));
	u = u+deltaU;
	x = Am*x+Bm*u+[0 1]'*d;
	buf = [buf; k u deltaU xe(3) x(2)];
	bufu = [bufu; [k k+1]' [u u]'];
end
k = buf(:, 1); y = buf(:, 5); u = buf(:, 2);
 ud = bufu(1:end-1, 2); ku = bufu(1:end-1, 1);
subplot 311, plot(k, y)
subplot 312, plot(k, u)
subplot 313, plot(ku, ud)