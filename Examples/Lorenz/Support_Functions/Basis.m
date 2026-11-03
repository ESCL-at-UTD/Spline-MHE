function B = Basis(x, knots, np)

x = x(:);

nBasis = length(knots) - np - 1;
nPts   = length(x);

% degree 0
N = zeros(nPts,nBasis+np);

for j = 1:(nBasis+np)

    N(:,j) = (x >= knots(j) & x < knots(j+1));

    % include right endpoint
    if knots(j+1) == knots(end)
        N(x==knots(end),j) = 1;
    end
end

% Cox-de Boor recursion
for k = 1:np

    Nk = zeros(nPts,nBasis+np-k);

    for j = 1:(nBasis+np-k)

        term1 = zeros(nPts,1);
        term2 = zeros(nPts,1);

        d1 = knots(j+k) - knots(j);
        if d1 > 0
            term1 = ((x-knots(j))/d1).*N(:,j);
        end

        d2 = knots(j+k+1) - knots(j+1);
        if d2 > 0
            term2 = ((knots(j+k+1)-x)/d2).*N(:,j+1);
        end

        Nk(:,j) = term1 + term2;
    end

    N = Nk;
end

B = N;
end