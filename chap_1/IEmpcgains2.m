function [F,phi,RRs] = IEmpcgains2(Ad,Bd,Cd,Np,Nc,rs)
%  for augmented system
RRs = rs*ones(Np,1);
phi = zeros(Np, Nc);
for a = 1:Np
    for b = 1:Nc
        if a < b
            phi(a,b) = 0;
        else
            phi(a,b) = Cd*Ad^(a-b)*Bd;
        end
    end
end
F = [];
for a = 1:Np
    F = [F; Cd*Ad^a];
end
