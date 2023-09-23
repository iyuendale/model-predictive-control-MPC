% Find th triplet matrices (A, B, C) in augmented system and eigenvalues

%  system matrices
Am = [1 1; 0 1]; Bm = [0.5 1]'; Cm = [1 0];

% triplet matrices of the augmented system
A = [Am zeros(size(Am, 1), size(Cm, 1));
	  Cm*Am eye(size(Cm, 1))]
B = [Bm; Cm*Bm]
C = [zeros(size(Cm, 1), size(Am, 1)) eye(size(Cm, 1))]

% eigenvalues of the system matrix
eigenvalue = eig(A)