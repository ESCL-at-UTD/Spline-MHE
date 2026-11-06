function dPhi = deltaBasis(t,Tau,np,delta)

if delta==0
    dPhi = Basis(t,Tau,np);
    return
end

Bm1 = deltaBasis(t,Tau,np-1,delta-1);

nBasis = length(Tau)-np-1;

dPhi = zeros(length(t),nBasis);

for j=1:nBasis

    d1 = Tau(j+np)-Tau(j);
    if d1~=0
        dPhi(:,j) = dPhi(:,j) + np/d1*Bm1(:,j);
    end

    d2 = Tau(j+np+1)-Tau(j+1);
    if d2~=0
        dPhi(:,j) = dPhi(:,j) - np/d2*Bm1(:,j+1);
    end

end
end