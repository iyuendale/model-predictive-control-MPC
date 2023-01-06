%% proobem from chapter 2
% problem 2.1
E= [3 1; 1 1]; F= [-15; -7];
% constraints
% -x1<=0, -x2<=0, x1+x2<=3, 2*x1- x2<=4, x2<=2
%% 1 = unconstrained minimization of J
% J= 1/2*(x'*E*x) + x'*F
% the general solution is: x = -E^(-1)*F
x_o = -inv(E)*F

%% 2 =  the last three constraints are violated

%% 3 = plot inequalities
v = -5:0.1:5;
[x1, x2]= meshgrid(v);
ineq = -x1 <= 0;
f = double(ineq);
subplot(1, 5, 1), surf(x1, x2, f); view(0, 90), axis square
xlabel 'x_1', ylabel ' x_2'
ineq2 = -x2 <= 0;
f2 = double(ineq2);
subplot(1, 5, 2), surf(x1, x2, f2), view(0, 90), axis square
xlabel 'x_1', ylabel ' x_2'
ineq3 = x1 + x2 <= 3; f3 = double(ineq3);
subplot(1, 5, 3), surf(x1, x2, f3), view(0, 90), axis square
xlabel 'x_1', ylabel ' x_2'
ineq4 = 2*x1 - x2 <= 4; f4 = double(ineq4);
subplot(1, 5, 4), surf(x1, x2, f4), view(0, 90), axis square
xlabel 'x_1', ylabel ' x_2'
ineq5 = x2 <= 2; f5 = double(ineq5);
subplot(1, 5, 5), surf(x1, x2, f5), view(0, 90), axis square
xlabel 'x_1', ylabel ' x_2'

% the first two constraints are active constraints as they are satisfied by
% the global solution, so taking the active constraints as
M_act = [-1 0; 0 -1]; b_act = [0 ; 0];
lambda_act = -(M_act*inv(E)*M_act')\(b_act + M_act*inv(E)*F)
x = -inv(E)*(F + M_act'*lambda_act)
% this gives an optimal solution that satisfies all of the constraints
