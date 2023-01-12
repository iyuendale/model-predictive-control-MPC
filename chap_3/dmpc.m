% Tuttorial 3.2
% function developed for generating the data matrices in the 
% design of DMPC where the cost function is 
% J = n'*E*n + 2*n'*H*x(ki)
% where n is vector of coefficients, E is omega and H is psi

function [E, H] = dmpc(A, B, a, N, Np, Q, R)
% A and B are found from the augmented system
% a contains the Laguerre pole locations(scaling factors) of each input
% N the number of terms forms for each input
% Np prediction horizon
% Q & R are weighting matrices on states and inputs
% R is assumed to be diagonal
% the cost function: J= eta ^T E eta +2 eta ^T H x(k_i)
[n, n_in] = size(B);
N_pa = sum(N);  % dimension of the eta(cofficient vector)
E = zeros(N_pa, N_pa); H = zeros(N_pa, n);
R_para = zeros(N_pa, N_pa);
n0 = 1; ne = N(1);
for i = 1:n_in-1
	R_para(n0:ne, n0:ne) = R(i, i)*eye(N(i), N(i));
	n0 = n0+N(i); ne = ne+N(i+1);
end
R_para(n0:N_pa, n0:N_pa) = R(n_in, n_in)*eye(N(n_in), N(n_in));
% initial conditions for convolution sum & calculate for i =1
S_in = zeros(n, N_pa);
[Al, L0] = lagd(a(1), N(1));
S_in(:, 1:N(1)) = B(:, 1)*L0';
In_s = 1;
for jj = 2:n_in
	[Al, L0] = lagd(a(jj), N(jj));
	In_s = N(jj-1)+In_s;
	In_e = In_s + N(jj) - 1;
	S_in(:, In_s:In_e) = B(:, jj)*L0';
end
S_sum = S_in;
phi = S_in;
E = (phi)'*Q*phi;
H = phi'*Q*A;
% iteration with respect to the prediction horizon
for i = 2:Np
	Eae = A^i;
	% calculate the finite sum S for each input
	%%%%
	% for each sample i
	%%%% calculate input number 1
	% specify the L0 and state matrix Al
	% associated with the first input
	[Al, L0] = lagd(a(1), N(1));
	% Laguerre function associated with input number 1
	S_sum(:, 1:N(1)) = A*S_sum(:, 1:N(1))+S_in(:, 1:N(1))*(Al^(i-1))';
	%%% move to input number 2 and so on
	In_s = 1;
	for kk = 2:n_in
		[Al, L0] = lagd(a(kk), N(kk));
		In_s = N(kk-1) + In_s;
		In_e = In_s+N(kk)-1;
		S_sum(:, In_s:In_e) = A*S_sum(:, In_s:In_e)+...
			S_in(:, In_s:In_e)*(Al^(i-1))';
	end
	phi = S_sum;
	E = E+phi'*Q*phi;
	H = H+phi'*Q*Eae;
end
E = E+R_para;

% This program generates the predictive control cost function for a system
% with an arbitrary number of inputs, states and outputs.