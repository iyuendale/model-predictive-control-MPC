% Tutorial 1.2
function [phi_phi,phi_F,phi_Rs,AA,BB,CC] = mpcgains(Am,Bm,Cm,Np,Nc,rs)
RRs = rs*ones(Np,1);
AA = [Am zeros(size(Am,1),size(Cm,1));Cm*Am eye(size(cd,1),size(cd,1))];
BB = [Bm; Cm*Bm]; CC = [zeros(size(Cm,1),size(Am,2)) eye(size(Cm,1))];

phi = [];
F = [];
Row = [];

for a = 1:Np
	F = [F; CC*AA^a];
	for b = 1:Nc
		if a < b
			Row = [Row zeros(size(CC, 1), size(BB, 2))];
		else
			Row = [Row CC*AA^(a-b)*BB];
		end
	end
	phi = [phi; Row]; Row = [];
end

phi_phi = phi'*phi;
phi_F = phi'*F;
phi_Rs = phi'*RRs;

