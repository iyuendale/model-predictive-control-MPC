Am= [1 .5 0; 0 1 -.1; 0 0 .8];
Bm= [.5 1 -.6]'; Cm = [1 0 1]; Bd = [1 0 0]';
Dm = 0;
%% Augmented system
A = [Am zeros(size(Am, 1), size(Cm, 1));
	  Cm*Am 1]
B = [Bm; Cm*Bm]
C = [zeros(1, size(Cm, 2)) 1]
D= 0;
eig(Am), eig(A)
%% Transfer function
% plant
sysp = ss(Am, Bm, Cm, Dm);
[num, den] = ss2tf(Am, Bm, Cm, Dm);
Gp = tf(num, den);
Zero_plant = zero(Gp)
Poles_plant= pole(Gp)
% augmented system
sysaug = ss(A, B, C, D);
[num, den] = ss2tf(A, B, C, D);
Gaug = tf(num, den);
Zeros_aug= zero(Gaug)
Poles_aug= pole(Gaug)
display('Transfer functions of plant and augmented system')
Gp, Gaug
