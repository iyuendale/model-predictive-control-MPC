A = [1 1 0; 0 1 0; 1 1 1]; B = [.5 1 .5]'; C = [0 0 1];
D = 0;
Nc = 5; Np = 30; rw = 10; rs = 1;

%% check controllability & observability
rank(ctrb(A, B)), rank(obsv(A, C))

%% MPC gains
[F,phi,RRs] = IEmpcgains2(A,B,C,Np,Nc,rs);
R = rw*eye(Nc, Nc);
Ky = inv(phi'*phi+R)*phi'*RRs;
Kmpc = inv(phi'*phi+R)*phi'*F;
% Kmpc = [0.8984 1.3521 0.4039]
%% observer
poles = [0.01 0.0105 0.011];
Kob = place(A', C', poles)'
%%
x = [0 0 0]'; xx = [0 0 .1]';
t = 0; h = 1; u = 0.4; buf = [t x' u u]; buf2 = [];
for k = 1:20
	deltaU = Ky - Kmpc*xx;
	u = u + deltaU(1); buf2 = [];
	for j = 1:Nc
		X = x; 
		X = A*X+B*deltaU(j);
		buf2 = [buf2; j+t X'];
		t2 = buf2(:, 1); x3 = buf2(:, 4);
		figure(2), plot(t2, x3, '--k'), hold on
% 		pause(.1)
	end
	xx = A*xx+B*deltaU(1)+Kob*C*(x-xx);
	x = A*x + B*deltaU(1);
	figure(2), plot(t+1, x(3), '*b', buf(:, 1), buf(:, 4), 'r')
	t = t+h; buf = [buf; t x' u deltaU(1)];
end
t = buf(:, 1); x1 = buf(:, 2); x2 = buf(:, 3);
x3 = buf(:, 4); u = buf(:, 5); deltaU = buf(:, 6);
figure (3)
subplot 311, plot(t, x3), legend 'y'
subplot 312, plot(t, u), legend 'u'
subplot 313, plot(t, deltaU), legend '\deltaU'

%% 
eig(A-B*Kmpc(1, :))
eig(A-Kob*C)