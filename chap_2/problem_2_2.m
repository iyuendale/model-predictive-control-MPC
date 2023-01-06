E = [3 1; 1 1]; F = [-15 -7]';
M = [-1 0; 0 -1; 1 1; 2 -1; 0 1];
b = [0 0 3 4 2]';
%% H and K
H = M*inv(E)*M'
K= b + M*inv(E)*F

%% global solution
x = -inv(E)*F
% fails to obey three of the constraints

%% Hildreth's algorithm
lambda = zeros(size(H, 1), 1);

for k = 1:50
	lambda_pre = lambda;
	for i = 1:size(lambda, 1)
		w = -(K(i) + H(i, :)*lambda - H(i, i)*lambda(i))/H(i, i);
		lambda(i) = max(0, w);
	end
	if abs(lambda-lambda_pre) < 0.000000001
		break;
	end
end

lambda
x = -inv(E)*(F + M'*lambda)