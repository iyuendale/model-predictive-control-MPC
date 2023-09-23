% Augmented design model

% continous-time system model
Ac = [0 1 0; 3 0 1; 0 1 0];
Bc = [1 1 3]'; Cc = [0 1 0]; Dc = 0;

% discretize the model
dt = 1;
[Ad, Bd, Cd, Dd] = c2dm(Ac, Bc, Cc, Dc, dt)

% form the augmented system
A = [Ad zeros(size(Ad, 1), size(Cd, 1));
	 Cd*Ad eye(size(Cd, 1))]
B = [Bd; Cd*Bd]
C = [zeros(size(Cd, 1), size(Ad, 1)) eye(size(Cd, 1))]