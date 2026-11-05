function [x_star,y_star,y_hat] = Lorenz_sim(true_param,N,dt,awg_n,plt)

%% Define a CasADi Model

import casadi.*

x0_star = true_param{1};
theta_star = true_param{2};

% Declare model variables
x = MX.sym('x', 3);
theta = MX.sym('theta', 3);

% Nonlinear Dynamics
xdot = [theta(1)*(x(2)-x(1))
        x(1)*(theta(2)-x(3)) - x(2)
        x(1)*x(2) - theta(3)*x(3)];

y = x(1);

fxdot = Function('fxdot', {x, theta}, {xdot});
fy = Function('fy', {x}, {y});

%% Formulate discrete time dynamics with fixed step Runge-Kutta 4 integrator
M = 4; % RK4 steps per interval
DT = dt/M;

X0 = MX.sym('X0', 3);
Theta = MX.sym('Theta', 3);
X = X0;
for j=1:M
   k1 = fxdot(X, Theta);
   k2 = fxdot(X + DT/2 * k1, Theta);
   k3 = fxdot(X + DT/2 * k2, Theta);
   k4 = fxdot(X + DT * k3, Theta);
   X = X + DT/6*(k1 +2*k2 +2*k3 +k4);
   Y = fy(X);
end
Fsim = Function('Fsim', {X0, Theta}, {X, Y},char('x0','theta'),char('x','y'));

%% Generate output trajectory
theta_k = DM(theta_star);
xk = DM(x0_star);
x_star = [xk];
y_star = [fy(x_star)];
for k=1:N
    ret = Fsim('x0',xk, 'theta',theta_k);
    xk = full(ret.x);
    yk = full(ret.y);
    x_star = [x_star xk];
    y_star = [y_star yk];
end

x_star = full(x_star);
y_star = full(y_star);


%% Corrupt output with Gaussian noise
rng(1);
y_hat = awgn(y_star,awg_n);


%% Plot generated Trajectory

if plt == 1
T = N*dt;
figure('Position', [50, 50, 1000, 350]);
tiledlayout(1,2,TileSpacing="compact");

nexttile
hold on
plot(linspace(0,T,N+1)',x_star(1,:),'-r',LineWidth=2)
plot(linspace(0,T,N+1)',x_star(2,:),'-b',LineWidth=2)
plot(linspace(0,T,N+1)',x_star(3,:),'-g',LineWidth=2)
xlim([0 T])
title("$\mathbf{States}$",Interpreter="latex")
xlabel('t(s)',Interpreter="latex")
ylabel('x',Interpreter="latex")
legend(["$x^{\star}_1$","$x^{\star}_2$","$x^{\star}_3$"],Interpreter="latex")
hold off
grid on

nexttile
hold on
plot(linspace(0,T,N+1)',y_hat,'k',LineWidth=2)
plot(linspace(0,T,N+1)',y_star,'r',LineWidth=2)
xlim([0 T])
title("$\mathbf{Output}$",Interpreter="latex")
xlabel('t(s)',Interpreter="latex")
ylabel('y',Interpreter="latex")
legend(["$\hat{\mathbf{Y}}$","$\mathbf{Y}^\star$"],Interpreter="latex")
hold off
grid on

end

end