% 
% E= eye(3, 3); F= [-2 -3 -1]';
% M= [1 1 1; 3 -2 -3; 1 -3 2]; b= [1 1 1]';
% eta= QPhild2(E, F, M, b)
function eta= QPhild(E, F, M, b)
% M & b= matrix and vector from expressions of constraints
[n1, ~]= size(M);
eta= -E\F;    % global solution
% kk= number of inactive constraints
kk = 0;
for i= 1:n1
		if (M(i, :)*eta > b(i))
			kk= kk+1;
		else
			kk= kk+0;
		end
end
if (kk == 0)
	return;    % if no inactive constraint is found
end

H= M*(E\M');
K= b+M*(E\F);
[n, m]= size(K);
lambda= zeros(n, m);       % initial

for km= 1:100
	lambda_pre= lambda;
	for i= 1:n
		w= -(K(i) + H(i, :)*lambda - H(i, i)*lambda(i,1))/H(i, i);
		lambda(i, 1) = max(0, w);
	end
	e= (lambda-lambda_pre)'*(lambda-lambda_pre);   % mean square error
	if (e < 10e-8)
		break;
	end
end
eta= eta- E\(M'*lambda);
