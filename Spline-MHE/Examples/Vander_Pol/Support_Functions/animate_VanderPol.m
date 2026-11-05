function animate_VanderPol(iteration_data,x_star,y_star,Y_hat,N,dt,nh,ns,n_x,n_y)

T = N*dt;
Ns = N/(ns*nh);

figure('Position', [50, 50, 1100, 300],'Color','white');

for i = 1:nh

iterData = iteration_data{i};

t_all = linspace(0,T/nh,(N/nh)+1)+((i-1)*(T/nh));
t = reshape(t_all(1:end-1),[],ns)';
tidx = find((rem(1:(N/nh),Ns)==0))+1;
t_spline = [t (t_all(tidx)'*NaN)];
t = reshape(t_spline',[],1)';

for j = 1:size(iterData,2)

%% State Plot

    subplot(1,3,1)
    hold on
    grid on
    
    xlim([0 T]);
    if i == 1
        lt1 = plot(linspace(0,T,N+1),x_star(1,:),"Color","#fc1303",LineWidth=2,LineStyle="--");
        lt2 = plot(linspace(0,T,N+1),x_star(2,:),"Color","#2803fc",LineWidth=2,LineStyle="--");
    end

    X_iter = [];
    for k = 1:ns
        X_iter = [X_iter [iterData{j}.x{k}(:,1:end-1) NaN*ones(n_x,1)]];
    end

    l1 = plot(t,X_iter(1,:),"Color","#fc1303",LineWidth=3,LineStyle="-");
    l2 = plot(t,X_iter(2,:),"Color","#2803fc",LineWidth=3,LineStyle="-");

    title("$\mathbf{States}$",Interpreter="latex")
    xlabel('t(s)',Interpreter="latex")
    ylabel('x',Interpreter="latex")
    legend([l1,l2],["$x_1$","$x_2$"],Interpreter="latex",NumColumns=3)
    hold off
    grid on

%% Output Plot

    subplot(1,3,2)
    hold on
    grid on
    
    xlim([0 T]);
    if i == 1
        lt4w = plot(linspace(0,T,N+1),Y_hat,"k",LineWidth=2,LineStyle="-");
        lt4 = plot(linspace(0,T,N+1),y_star,"Color","#419100",LineWidth=2,LineStyle="--");
    end
    
    Y_iter = [];
    for k = 1:ns
        Y_iter = [Y_iter [iterData{j}.y{k}(:,1:end-1) NaN*ones(n_y,1)]];
    end

    l4 = plot(t,Y_iter,"Color","#419100",LineWidth=2,LineStyle="-");

    title("$\mathbf{Output}$",Interpreter="latex")
    xlabel('t(s)',Interpreter="latex")
    ylabel('y',Interpreter="latex")
    legend([lt4w,l4],["$\mathbf{\hat{Y}}$" "$\mathbf{Y}$"],Interpreter="latex",NumColumns=2)
    hold off
    grid on

%% State-Space Plot

    subplot(1,3,3)
    hold on
    grid on
    if i == 1
        lt_sp = plot(x_star(1,:),x_star(2,:),'k-',LineWidth=2,LineStyle="--");
    end
    l_sp = plot(X_iter(1,:),X_iter(2,:),'r-',LineWidth=2,LineStyle="-");

    title("$\mathbf{State-Space}$",Interpreter="latex")
    xlabel('$x_1$',Interpreter="latex")
    ylabel('$x_2$',Interpreter="latex")
    hold off
    grid on

%% Animate

    pause(0.01)
    drawnow

    if j ~= size(iterData,2)
        delete(l1);
        delete(l2);
        delete(l4);
        delete(l_sp);
    else
        L1 = l1;
        L2 = l2;
        L4 = l4;
        Lsp = l_sp;
    end

end


end



end

