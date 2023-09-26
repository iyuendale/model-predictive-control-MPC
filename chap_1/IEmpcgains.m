% Tutorial 1.2
function [F,E,RRs,AA,BB,CC] = IEmpcgains(Ad,Bd,Cd,Np,Nc,rs)
RRs = rs*ones(Np,1);
AA = [Ad zeros(size(Ad,1),size(Cd,1));Cd*Ad eye(size(cd,1),size(cd,1))];
BB = [Bd; Cd*Bd]; CC = [zeros(size(Cd,1),size(Ad,2)) eye(size(Cd,1))];
F = zeros(Np, Nc);
for a = 1:Np
    for b = 1:Nc
        if a < b
            F(a,b) = 0;
        else
            F(a,b) = CC*AA^(a-b)*BB;
        end
    end
end
E = [];
for a = 1:Np
    E = [E; CC*AA^a];
end
