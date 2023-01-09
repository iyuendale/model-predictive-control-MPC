% creates the Laguerre network matrix A_l and initial vector L0
function [A_l, L0] = lagd(a, N)
v(1, 1) = a; L0(1, 1) = 1;
for k = 2:N
	v(k, 1) = (-a).^(k-2)*(1-a*a);
	L0(k, 1) = (-a).^(k-1);
end

L0 = sqrt((1-a*a))*L0;
A_l(:, 1) = v;
for i = 2:N
	A_l(:, i) = [zeros(i-1, 1); v(1: N-i+1, 1)];
end
