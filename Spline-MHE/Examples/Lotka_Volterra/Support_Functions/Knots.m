function Tau = Knots(np,n_tau,t)

Tau = linspace(min(t(1,:)), max(t(1,:)), n_tau);
Tau = augknt(Tau, np + 1);

end