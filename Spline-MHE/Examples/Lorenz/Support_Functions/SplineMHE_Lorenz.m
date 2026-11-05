% function Sol = SplineMHE_Lorenz(C_init,Cest_init,xi0,Y_obs_all,Y_true,spline_param,N,n_iter,n_mhe,dt,Wx,Wy)
function Sol = SplineMHE_Lorenz(x0_init,theta_init,Y_hat_est,delta,Kappa,N,dt,ns,nh,np,n_tau)

global iteration_data Sol_iter Sol;
iteration_data = cell(1,nh);


Ns = N/(ns*nh); % number of discrete points within a spline segment
N = N/nh; % number of discrete points within a horizon

n_x = 3; % number of states
n_y = 1; % number of outputs
n_theta = 3; % number of parameters


%% Construct spline basis matrices
tidx = find((rem(1:N,Ns)==0))+1;
idx = [reshape(1:ns*(Ns),[],ns)' tidx'];

t_all = linspace(0,N,N+1);
t = reshape(t_all(1:end-1),[],ns)';
t = [t t_all(tidx)'];

Tau = Knots(np,n_tau,t(1,:)); % knots  for spline segments

Phi_hat = Basis(t(1,:), Tau, np); % B-spline basis for spline segments
dPhi_hat = {};
for d=1:delta
    dPhi_hat{end+1} = deltaBasis(t(1,:),Tau,np,d); % B-spline derivative basis for spline segments
end

Phi_hat_dagger0 = (Phi_hat'*Phi_hat)\Phi_hat';



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Spline-MHE NLP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

opti = casadi.Opti(); % Optimization problem

%% NLP decision variable and parameter definition

C_bar = {};
for i = 1:ns
    C_bar{end+1} = opti.variable((np+(n_tau-1)),n_x); % state trajectory smoothing spline scaling coefficients
end

C = {};
for i = 1:ns
    C{end+1} = opti.variable((np+(n_tau-1)),n_y); % output trajectory smoothing spline scaling coefficients
end

theta = opti.variable(n_theta,1); % parameter vector

% Y_est = opti.parameter(n_y,N+1); % observed output smoothing spline estimate

C_hat = {};
for i = 1:ns
    C_hat{end+1} = opti.parameter((np+(n_tau-1)),n_y); % output trajectory smoothing spline scaling coefficients
end


%% Lorenz system nonlinear dynamics

fxdot = @(x,theta)[theta(1).*(x(2,:)-x(1,:)); x(1,:).*(theta(2)-x(3,:)) - x(2,:); x(1,:).*x(2,:) - theta(3)*x(3,:)];
fy = @(x) [x(1,:)];


%% Spline-MHE Iterations

X_cost = [];
Y_cost = [];
dY_cost = [];

for i = 1:ns

    x = (Phi_hat*C_bar{i})';
    dx = (dPhi_hat{1}*C_bar{i})';

    % xdot constraints
    dX = fxdot(x,theta);
    opti.subject_to(dx(:,Kappa)./dt == dX(:,Kappa));

    % output constraints
    Y = fy(x);
    y = (Phi_hat*C{i})';
    opti.subject_to(y(:,Kappa) == Y(:,Kappa));

    % Loss function
    Y_cost = [Y_cost sum((Phi_hat*C{i}-Phi_hat*C_hat{i})'*(Phi_hat*C{i}-Phi_hat*C_hat{i}))];
    for d = 1:delta
        dY_cost = [dY_cost sum((dPhi_hat{d}*C{i}-dPhi_hat{d}*C_hat{i})'*(dPhi_hat{d}*C{i}-dPhi_hat{d}*C_hat{i}))];
    end

    % Spline segment endpoint penalization
    if i ~= ns
        x_end = (Phi_hat*C_bar{i+1})';
        X_cost = [X_cost sum((x(:,end) - x_end(:,1)).^2)];
    end
end


obj = sum(Y_cost) + sum(dY_cost) + sum(X_cost);
opti.minimize(obj);

t_construct = toc;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solve Spline-MHE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sol = {};
Sol_iter = {};

% Lorenz system simulation at initial guess
theta_min = theta_init; % Initializing theta_min at initial parameter guess
init_param = {x0_init,theta_min}; % Initializing states and parameters
[x_sim,y_sim,~] = Lorenz_sim(init_param,N,dt,0,false);

%% Moving Horizon Iterations
for i = 1:nh
    
    % Initialize C_bar, C and C_hat for ns spline segments
    C_bar_init = {};
    for j = 1:ns
        C_bar_init{end+1} = Phi_hat_dagger0*x_sim(:,idx(j,:))'; % simulation state trajectory smoothing spline segment scaling coefficient estimate
    end
    
    C_init = {};
    for j = 1:ns
        C_init{end+1} = Phi_hat_dagger0*y_sim(:,idx(j,:))'; % simulation output trajectory smoothing spline segment scaling coefficient estimate
    end
    
    C_hat_init = {};
    for j = 1:ns
        C_hat_init{end+1} = Phi_hat_dagger0*Y_hat_est(:,idx(j,:)+((i-1)*N))'; % measured output trajectory smoothing spline segment scaling coefficient estimate
    end

    % Assigning initial values to NLP variables and parameters
    for j = 1:ns
        opti.set_initial(C_bar{j}, C_bar_init{j});
        opti.set_initial(C{j}, C_init{j});
        opti.set_value(C_hat{j}, C_hat_init{j});
    end
    opti.set_initial(theta, theta_min);

    % Solve Spline-MHE NLP
    opti.callback(@(l) append_data(l, opti.debug.value(opti.x),n_x,n_y,n_theta,ns,np,n_tau,Phi_hat));
    ipopt_options = struct('print_level', 1, 'tol', 1e-5);
    solver_opts = struct('ipopt', ipopt_options); 
    opti.solver('ipopt',solver_opts);
    sol = opti.solve();
    
    % Parse solution
    x_sol = [];
    for j = 1:ns
        x_sol = [x_sol (Phi_hat*sol.value(C_bar{j}))'];
    end
    
    y_sol = [];
    for j = 1:ns
        y_sol = [y_sol (Phi_hat*sol.value(C{j}))'];
    end

    theta_sol = sol.value(theta);

    x0_sol = x_sol(:,1);
    Cost = sol.value(obj)/(N+1);
    stats = sol.stats();
    
    Sol{end+1} = {x0_sol;theta_sol;Cost;stats.iter_count;t_construct;stats.t_wall_total};
    iteration_data{i} = Sol_iter;
    Sol_iter = {};

    % Reset Moving Horizon Iteration
    if i ~= nh

        x0_sol = x_sol(:,end);
        % Lorenz system simulation at theta_min
        init_param = {x0_sol,theta_min}; % update states and parameters for next horizon
        [x_sim_min,y_sim_min,~] = Lorenz_sim(init_param,N,dt,0,false);

        % Lorenz system simulation at new solution
        init_param = {x0_sol,theta_sol}; % update states and parameters for next horizon
        [x_sim_sol,y_sim_sol,~] = Lorenz_sim(init_param,N,dt,0,false);

        % Convergence costs
        J_theta_min = (1/(Ns+1))*(sum((y_sim_min-Y_hat_est(:,(1:N+1)+(i*N)))'*(y_sim_min-Y_hat_est(:,(1:N+1)+(i*N))))); % cost at theta_min
        J_theta_sol = (1/(Ns+1))*(sum((y_sim_sol-Y_hat_est(:,(1:N+1)+(i*N)))'*(y_sim_sol-Y_hat_est(:,(1:N+1)+(i*N))))); % cost at theta_sol

        % Update theta_min, x_sim, and y_sim
        if J_theta_sol <= J_theta_min
            theta_min = theta_sol;
            x_sim = x_sim_sol;
            y_sim = y_sim_sol;
        else
            theta_min = theta_min;
            x_sim = x_sim_min;
            y_sim = y_sim_min;
        end
    end
end

Sol = cell2table([Sol{:,:}],RowNames=["x*(0)" "θ*" "Convergence Cost" "Solver Iterations" "NLP Construction Time (s)" "NLP Solve Time (s)"],VariableNames=["Horizon-"+string(1:nh)]);

end



%% Support function

function append_data(iteration, variable, n_x, n_y, n_theta, ns, np, n_tau, Phi_hat)
    global Sol_iter;
    sol = full(variable);

    C_bar = sol(1:((np+(n_tau-1))*n_x*ns));
    C = sol(((np+(n_tau-1))*n_x*ns)+1:end-n_theta);

    C_bar = reshape(C_bar,np+(n_tau-1),[]);
    C = reshape(C,np+(n_tau-1),[]);

    x_sol = {};
    j = 1;
    for i = 1:ns
        x_sol{end+1} = (Phi_hat*C_bar(:,j:j+(n_x-1)))';
        j = j+n_x;
    end

    y_sol = {};
    j = 1;
    for i = 1:ns
        y_sol{end+1} = (Phi_hat*C(:,j:j+(n_y-1)))';
        j = j+n_y;
    end

    theta = sol(end-(n_theta-1):end);

    Sol_iter{end+1} = struct('iter', iteration, 'x', {x_sol}, 'y', {y_sol}, 'theta', {theta});
end