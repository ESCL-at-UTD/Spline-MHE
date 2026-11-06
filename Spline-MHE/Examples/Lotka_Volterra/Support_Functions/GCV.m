function [lambda,Y_hat_est] = GCV(Y_hat,delta,np,n_tau,N,dt,plt)

t = linspace(0,N,N+1); % discrete time steps 
Tau = Knots(np,n_tau,t); % knots

Phi = Basis(t, Tau, np); % B-spline basis
dPhi = deltaBasis(t,Tau,np,delta); % B-spline derivative basis


lambda = 1; % inital lambda

Pseudo_Phi = (Phi'*Phi + lambda*((dPhi'./(dt^delta))*(dPhi./(dt^delta))))\Phi';

% argmin estimate of lambda

GCV_est1=1;
GCV_est2=10;

while abs(GCV_est1-GCV_est2)>=1e-3 % Recursive estimate for GCV minimizing lambda

    lambda1 = lambda/2;
    lambda2 = lambda*2;
    
    GCV_est1 = GCVest(Phi,dPhi,Y_hat,lambda1,delta,N,dt);
    GCV_est2 = GCVest(Phi,dPhi,Y_hat,lambda2,delta,N,dt);

    if GCV_est1<=GCV_est2
        lambda = lambda1;
    else
        lambda = lambda2;
    end
end

% Smoothing spline estimate of output data
Pseudo_Phi = (Phi'*Phi + lambda*((dPhi'./(dt^delta))*(dPhi./(dt^delta))))\Phi';
Y_hat_est = (Phi*Pseudo_Phi*Y_hat')';


% Plot estimate

if plt == true
figure()
T = N*dt;
hold on
plot(linspace(0,T,N+1)',Y_hat,'k',LineWidth=2)
plot(linspace(0,T,N+1)',Y_hat_est,'r',LineWidth=2)
xlim([0 T])
title("$\mathbf{Smoothing \ Spline \ Estimate}$",Interpreter="latex")
xlabel('t(s)',Interpreter="latex")
ylabel('y',Interpreter="latex")
legend(["$\hat{\mathbf{Y}}$","$\hat{\mathbf{Y}}^{est}$"],Interpreter="latex")
hold off
grid on
end

end


%% GCV Function
function GCV_est = GCVest(Phi,dPhi,Y_hat,lambda,delta,N,dt)

Pseudo_Phi = (Phi'*Phi + lambda*((dPhi'./(dt^delta))*(dPhi./(dt^delta))))\Phi';
Y_hat_est = (Phi*Pseudo_Phi*Y_hat')';

GCV_est = trace((Y_hat_est-Y_hat)'*(Y_hat_est-Y_hat))/((N+1)*(1-(1/(N+1))*trace(Phi*Pseudo_Phi))^2);

end