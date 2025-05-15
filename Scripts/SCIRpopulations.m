function SCIRpopulations()
    param.numvar = 4; % do not change

    % initial parameters for SIR
    param.P     = 4;                                    % number of populations
    % to change P, also change the size of N, beta, and gamma
    param.N     = [100;100;100;100];                    % initial population sizes
    param.beta  = 2.0*[1.00;1.00;1.00;1.00]./param.N;   % infection rates
    param.gamma = 1.0*[1.0;1.0;1.0;1.0];                % recovery rates
    param.nu    = 0.01;                                 % base travel rate
    
    % add parameters for SCIR
    param.alpha = 1.0*[0.25;0.25;0.25;0.25]./param.N;   % transmission rate from carriers to susceptibles
    param.alphapr = 0./param.N;                         % reinfection rate of recovered individuals by carriers
    param.betapr = 0./param.N;                          % reinfection rate of recovered individuals by infected individuals
    param.omega = 0.5;                                  % recovery rate of carriers
    param.delta = 0.3;                                  % transition rate from carrier to infected
    
    % setup the travel (immigration / emigration) rates for each compartment
    % T(i,j) = rate from population j to i
    % the diagonal of T is effectively ignored
    % here are some simple network structures

    % no travel
    Tnone  = zeros(param.P,param.P);
    % travel between all populations
    Tall   = ones(param.P,param.P);
    % travel on a chain: 1 <-> 2 <-> 3 <-> ... <-> P
    Tchain = diag(ones(param.P-1,1),1)+diag(ones(param.P-1,1),-1);
    % travel only to and from population 1
    Tonly1 = [0 ones(1,param.P-1); ones(param.P-1,1) zeros(param.P-1,param.P-1)];
    % travel in one direction in a circle: 1 <- 2 <- 3 <- ... <- P <- 1
    Tcir   = diag(ones(param.P-1,1),1);
    Tcir(param.P,1) = 1;
    % Watts-Strogatz small world network
    %comment out the next line if get an error "undefined function repelem"
    Tsw = WattsStrogatz(param.P,1,0.5); % this is randomly generated each time
    
    param.T_S = param.nu*Tall;
    param.T_I = param.nu*Tnone;
    param.T_R = param.nu*Tall;
    param.T_C = param.nu*Tall;
    
    % visualize the networks
    shownetwork = true; % set to false if you get an error with "digraph"
    if shownetwork
        lab = {'S','I','R', 'C'};
        ids = cell(param.P*param.numvar,1);
        for c=1:param.numvar
            for i=1:param.P
                ids{(c-1)*param.P+i} = [ num2str(i) ' (' lab{c} ')'];
            end
        end
        % travel matrix
        M = blkdiag(param.T_S, param.T_I, param.T_R, param.T_C);
        M(1:size(M,1)+1:end) = 0;
        bg = digraph(M,ids);
        figure();
        plot(bg)
    end
    
    % simulation parameters
    tf = 20;
    numtimepts = 500;
    showpieannimation = false;
    
    % initial condition
    % everyone is susceptible, except a small percentage in the last population
    SIR0 = [param.N zeros(param.P,1) zeros(param.P,1) zeros(param.P,1)];
    SIR0(end,:) = param.N(end) * [0.95, 0.05, 0.00, 0.00];

    % compute the (initial) basic reproductive number for each population 
    R0 = (param.alpha .* (param.gamma + param.delta) + param.beta .* (param.omega + param.delta)) ...
        .* SIR0(:,1) ./ (param.gamma .* (param.omega + param.delta));
    
    % simulate
    [t,SIR] = ode45( @(t,x) SIR_rhs(t,x,param), linspace(0,tf,numtimepts), SIR0(:) );
    SIR = reshape(SIR,[numel(t),param.P,param.numvar] );
    S = SIR(:,:,1); I = SIR(:,:,2); R = SIR(:,:,3); C = SIR(:,:,4);
    
    % check for conservation
    if any(abs( sum(S+I+R+C,2) - sum(param.N) ) > 1e-8)
        error('Numbers not conserved over all populations')
    end
    
    % plot
    fh = figure();
    set(fh,'Position',[300 100 600 800]);
    clrs = [0 0 240; 220 20 60; 60 60 60; 100 100 100] / 256;
    clab = {'S','I','R','C'};
    for i=1:param.P
        for c=1:param.numvar
            labs{(i-1)*param.numvar+c} = [clab{c} '_' num2str(i)];
        end
    end
    
    % plot the time courses as area plots
    subplot(4,1,1); hold on;
    area(t,reshape(permute(SIR,[1,3,2]),numtimepts,param.P*param.numvar))
    xlabel('time')
    ylabel('population size')
    axis([0 tf 0 sum(param.N)])
    title(['(inital) R_0 = ' sprintf('%0.1f  ',R0')])
    colororder(clrs)
    
    % plot the pie chart distribution at each time
    tol = 1e-3*sum(param.N);
    subplot(4,1,2:4)
    cm = repmat(clrs,param.P,1);
    if showpieannimation
        tivals = 1:numtimepts;
    else
        tivals = numtimepts;
    end
    for ti=tivals
        figure(fh); cla;
        SIRcur = [S(ti,:);I(ti,:);R(ti,:);C(ti,:)];
        SIRcur = SIRcur(:);
        bigenough = SIRcur>tol;
        SIRcur = SIRcur(bigenough);
        ph = pie(SIRcur(:),labs(bigenough));
        ind = find(bigenough);
        for j=1:2:numel(ph)
            ph(j).FaceColor = cm(ind((j+1)/2),:);
        end
        title(['t = ' num2str(t(ti),'%0.1f')])
        drawnow
    end
    
    % display the final global distribution
    disp('(% of global population)')
    disp('i    S     I     R     C')
    for i=1:param.P
        fprintf('%d', i);
        for v=1:param.numvar
            fprintf(' % 5.1f', SIR(end,i,v)/sum(param.N)*100);
        end
        fprintf('\n');
    end
    disp('   ----  ----  ----')
    fprintf('  ')
    fprintf('% 5.1f ',squeeze(sum(SIR(end,:,:),2))'/sum(param.N)*100)
    fprintf('\n')

end


function dxdt = SIR_rhs(t,x,param)
    % Reshape x back to the matrix where each column represents a compartment
    X = reshape(x, [], param.P);

    % Extract compartment values for each population
    S = X(:,1); % Susceptible
    I = X(:,2); % Infected
    R = X(:,3); % Recovered
    C = X(:,4); % Carrier

    % Compute the vaccination rate depending on time t for each population
    % Ensuring that vaccination does not exceed the number of susceptible individuals
    v_t = 0; %min(100 * log(t + 1) ./ param.N, S);

    % Calculate derivatives
    dSdt = -param.beta .* S .* (I + C) - param.alpha .* S .* (I + C) - v_t + param.T_S*S - sum(param.T_S,1)'.*S;
    dIdt = param.beta .* S .* (I + C) + param.betapr .* R .* (I + C) - param.gamma .* I + param.delta .* C + param.T_I*I - sum(param.T_I,1)'.*I;
    dRdt = param.gamma .* I + param.omega .* C - param.alphapr .* R .* (I + C) - param.betapr .* R .* (I + C) + v_t + param.T_R*R - sum(param.T_R,1)'.*R;
    dCdt = param.alpha .* S .* (I + C) + param.alphapr .* R .* (I + C) - param.omega .* C - param.delta .* C + param.T_C*C - sum(param.T_C,1)'.*C;

    dxdt = [dSdt; dIdt; dRdt; dCdt];
end

