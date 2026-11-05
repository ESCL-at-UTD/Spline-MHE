function Phi = Basis(t, Tau, np)

t = t(:);

nBasis = length(Tau) - np - 1;
nPts   = length(t);

% degree 0
N = zeros(nPts,nBasis+np);

for j = 1:(nBasis+np)

    N(:,j) = (t >= Tau(j) & t < Tau(j+1));

    % include right endpoint
    if Tau(j+1) == Tau(end)
        N(t==Tau(end),j) = 1;
    end
end

% Cox-de Boor recursion
for k = 1:np

    Nk = zeros(nPts,nBasis+np-k);

    for j = 1:(nBasis+np-k)

        term1 = zeros(nPts,1);
        term2 = zeros(nPts,1);

        d1 = Tau(j+k) - Tau(j);
        if d1 > 0
            term1 = ((t-Tau(j))/d1).*N(:,j);
        end

        d2 = Tau(j+k+1) - Tau(j+1);
        if d2 > 0
            term2 = ((Tau(j+k+1)-t)/d2).*N(:,j+1);
        end

        Nk(:,j) = term1 + term2;
    end

    N = Nk;
end

Phi = N;
end