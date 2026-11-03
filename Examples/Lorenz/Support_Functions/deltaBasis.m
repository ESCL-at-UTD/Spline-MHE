function D = deltaBasis(x,knots,np,delta)

if delta==0
    D = Basis(x,knots,np);
    return
end

Bm1 = deltaBasis(x,knots,np-1,delta-1);

nBasis = length(knots)-np-1;

D = zeros(length(x),nBasis);

for j=1:nBasis

    d1 = knots(j+np)-knots(j);
    if d1~=0
        D(:,j) = D(:,j) + np/d1*Bm1(:,j);
    end

    d2 = knots(j+np+1)-knots(j+1);
    if d2~=0
        D(:,j) = D(:,j) - np/d2*Bm1(:,j+1);
    end

end
end