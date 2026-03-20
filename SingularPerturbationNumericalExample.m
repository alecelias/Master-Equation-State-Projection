clear all
% close all
clc

global n W S Q QT FNorm invFNorm indxpiv epsl

%% Set up the problem
% Define the number of states and the propensities for the allowed transitions
% kfij refers to transition from state i to j ( i => j )

RunParSet1
% RunParSet2
% RunParSet3
disp(['Number of states: ' num2str(m)]);

% epsl_range = [3.2 1.6 0.8 0.4 0.2 0.1 0.05 0.025 0.0125 0.00625 0.003125];
% epsl_range = logspace(-0.1,-3.1,30);
epsl_range = [0.5];
normErr0 = [];
normErr1 = [];
normErr2 = [];
maxNormalizErrFullSys = [];
maxNormalizErrReduSys0 = [];
maxNormalizErrReduSys1 = [];
maxNormalizErrReduSys2 = [];

%% Flags for plotting
% plotVerbose creates plots for each individual approximation as well as
% the full solution.
% fariplot creates plots showing the solution to the full equation as well
% as the different approximate solutions.
% fullplot creates plots that show the full solution for each of the
% individual states, as well as the superbasins that they belong to.
% approxplot creates plots that show the magnitude of each term in the
% perturbation expansion.

plotVerbose = false;
% plotVerbose = true;

fairplot = false;
% fairplot = true;

fullplot = false;
% fullplot = true;

approxplot = false;
% approxplot = true;

%% Set up the transition matrix

F = zeros(m,m);
S = zeros(m,m);
for i = 1:m
    for j = 1:m
        kstring = ['kf' num2str(i) num2str(j)];
        if exist(kstring,'var')
            % disp(['Adding contributions from ' kstring ' = ' num2str(eval(kstring))])
            % Probability gain term for state j <= i
            F(j,i) = F(j,i) + eval(kstring);
            % Probability loss term for state i => j
            F(i,i) = F(i,i) - eval(kstring);
            % disp(F)
        end
        kstring = ['ks' num2str(i) num2str(j)];
        if exist(kstring,'var')
            % disp(['Adding contributions from ' kstring ' = ' num2str(eval(kstring))])
            % Probability gain term for state j <= i
            S(j,i) = S(j,i) + eval(kstring);
            % Probability loss term for state i => j
            S(i,i) = S(i,i) - eval(kstring);
            % disp(S)
        end

    end
end

%% Create necessary matrices

Nsp = ceil(null(transpose(F))); % left nullspace of F
n = rank(Nsp); % number of superbasins
QT = Nsp;
% QT = QT(:,n:-1:1); % rearrange columns (for convenient relabelling of superbasins)
% QT = QT(:,[3 1 4 2]); % some arbitrary permutation of columns
Q = transpose(QT);
QTQ = QT*Q;
[REFQTQ,piv] = rref(QTQ); % row echelon form of QT*Q
indxpiv = false(1,m);
indxpiv(piv) = true;
% Create a full-rank matrix F with the normalisation equations for each
% superbasin
FNorm = F;
FNorm(indxpiv,:) = QTQ(indxpiv,:);
invFNorm = inv(FNorm);

%% Solve the full and reduced models
% Set options for the solver and the time and timesteps.

% ODESolOptions = odeset('RelTol',1e-7,'AbsTol',1e-10); % Options for the ODE solver
ODESolOptions = odeset('RelTol',3e-14,'AbsTol',3e-14); % Options for the ODE solver
tspan = [0 10];
t = linspace(tspan(1),tspan(2),1001); % times to evaluate the solution
Dt = t(2)-t(1);

%% Create initial conditions

% p = [0.3 0.5 0.2].';
p = [0.1 0.05 0.15 0.07 0.03 0.25 0.15 0.2].';
rhstmp = QT*Q*p;
rhstmp(~indxpiv) = 0;
pini0 = invFNorm*rhstmp; % use quasi equilibrated (within superbasins) probabilities

ptilde012ini = [Q*pini0; zeros(2*n,1)];

%% Solve the approximate system
% The solution to this does not change with epsilon so this can be done
% outside of the loop

solRedu012 = ode15s(@ReduOe2SystemRHS,tspan,ptilde012ini,ODESolOptions);
psolRedu012 = deval(solRedu012,t);

for epsl = epsl_range

    W = 1/epsl*F + S;
    % disp(W)

    % %% Cross-check our formulas
    % pini = zeros(m,1);
    % pini(1) = 1/3;
    % pini(2) = 2/3;
    % 
    % % Solve for a few initial timesteps with forward Euler
    % p = pini;
    % Dt = 0.00005;
    % for i = 1:1000
    %     p = p + W*p*Dt;
    %     % stem(p)
    %     % title(['step ' num2str(i*Dt)])
    %     % drawnow
    % end
    % if any(isnan(p))
    %     disp(p)
    %     error('At least one of the elements of p has become NaN!')
    % end

    % %% Demonstrate the solution of the O(1/eps) equations
    % rhstmp = QTQ*p; % right-hand side must have the values of ptilde
    % % in the rows corresponding to the normalisation equations
    % rhstmp(~indxpiv) = 0; % the other rows should have zeros (fluxes between superbasin's states balance out)
    % pQEq = invFNorm*rhstmp; % this gives the probabilities quasi-equilibrated within each superbasin
    % sum(pQEq);
    % % Validate the projector: the following expression should evaluate to the
    % % zero matrix
    % Q*F;
    

    % %% Check the validity of the distribution matrix K
    % K = diag(p)*inv(diag(QT*Q*p))*QT;
    % ptilde = Q*p;
    % p - K*ptilde; % this should give the zero vector (up to numerical accuracy)

    %% Create initial vector that is quasi-equilibrated.

    dpinitilde0dt = ReduOe0SystemRHS(0,Q*pini0);
    
    rhstmp = QT*dpinitilde0dt;
    rhstmp(~indxpiv) = 0;
    dpini0dt = invFNorm*rhstmp;
    
    rhstmp = QT*zeros(n,1);
    rhstmp2 = dpini0dt - S*pini0;
    rhstmp(~indxpiv) = rhstmp2(~indxpiv);
    pini1 = invFNorm*rhstmp;
    
    dpinitilde1dt = ReduOe1SystemRHS(0,[Q*pini0; Q*pini1]);
    dpinitilde1dt = dpinitilde1dt(n+1:2*n);
    
    rhstmp = QT*dpinitilde1dt;
    rhstmp(~indxpiv) = 0;
    dpini1dt = invFNorm*rhstmp;
    
    rhstmp = QT*zeros(n,1);
    rhstmp2 = dpini1dt - S*pini1;
    rhstmp(~indxpiv) = rhstmp2(~indxpiv);
    pini2 = invFNorm*rhstmp;
    
    pini = pini0 + epsl*pini1 + epsl^2*pini2;
    % pini = pini0;
    % pini = p;
    % pini = pini0 + epsl * pini1;

    %% Solve the full system
    solFull = ode15s(@FullSystemRHS,tspan,pini,ODESolOptions);
    psolFull = deval(solFull,t);
    maxNormalizErrFullSys = [maxNormalizErrFullSys max(sum(psolFull,1)-1)];

    %% Calculate approximate solutions for given epsilon and evaluate error
    psolRedu2 = psolRedu012(1:n,:) + epsl*psolRedu012(n+1:2*n,:) ...
        + epsl^2*psolRedu012(2*n+1:3*n,:);
    psolRedu1 = psolRedu012(1:n,:) + epsl*psolRedu012(n+1:2*n,:);
    psolRedu0 = psolRedu012(1:n,:);

    maxNormalizErrReduSys0 = [maxNormalizErrReduSys0 max(sum(psolRedu0,1)-1)];
    maxNormalizErrReduSys1 = [maxNormalizErrReduSys1 max(sum(psolRedu1,1)-1)];
    maxNormalizErrReduSys2 = [maxNormalizErrReduSys2 max(sum(psolRedu2,1)-1)];

    %% Evaluate the error of all the systems
    Err0 = sqrt(sum((Q*psolFull-psolRedu0).^2*Dt,2));
    normErr0 = [normErr0 norm(Err0)];

    Err1 = sqrt(sum((Q*psolFull-psolRedu1).^2*Dt,2));
    normErr1 = [normErr1 norm(Err1)];

    Err2 = sqrt(sum((Q*psolFull-psolRedu2).^2*Dt,2));
    normErr2 = [normErr2 norm(Err2)];
    
    %% Plotting if needed

    if fullplot
        fpV1 = figure;
        % fpV1.Position
        % fpV1.Position = [1, 1, 960, 720];
        % fontsize(fpV1, 20, "points")
        
        hold on
        b = plot(t, Q*psolFull,'--', LineWidth=1);
        C = colororder;
        D = [C(1:4, :); C(3, :); C(3, :); C(3, :); C(2, :); C(2, :); C(1, :); C(1, :); C(4, :)];
        colororder([D]);
        a = plot(t,psolFull, '-', LineWidth=2);
        % set(gca,'xscale','log')
        % title('Full Model - Full Solution')
        % legend('1', '2', '3', '4', '5', '6', '7', '8');
        xlabel('t');
        ylabel('Probability');
        title('\epsilon = 0.5')
        ylim([0, 1]);
        hleglines = [a(1), a(4), a(6), a(8)];
        legend(hleglines, 'Superbasin 1', 'Superbasin 2', 'Superbasin 3', 'Superbasin 4')
        fontsize(fpV1, scale=1.5);
        box on
        hold off
        exportgraphics(fpV1, 'Full-solution(8)05.pdf')
    end

    if approxplot
        f1 = figure;
        a = plot(t, psolRedu012(1:n, :), '-', LineWidth=2);
        C = colororder;
        colororder(C(1:4, :));
        hold on
        b = plot(t, psolRedu012(n+1:2*n, :), '--', LineWidth=2); 
        c = plot(t, psolRedu012(2*n+1:3*n, :), ':', LineWidth=2);
        colororder({'black', 'black'});
        yyaxis left
        ylabel('Probability');
        yyaxis right
        ylabel('Probability/\epsilon, Probability/\epsilon^2');
        yline(0, '--');
        hleglines = [a(3), b(3), c(3)];
        lll = legend(hleglines, '$\tilde p^{(0)}$', '$\tilde p^{(1)}$', '$\tilde p^{(2)}$', 'Interpreter', 'latex', 'Location', 'east');
        ylim([-0.2, 0.8]);
        xlabel('t');
        hold off
        % ylabel('Probability');
        fontsize(f1, scale=1.5);
        lgps = lll.Position
        lll.Position = lgps + [0, 0.05, 0, 0]
        exportgraphics(f1, '8stateexpansionterms.pdf');

        f2 = figure;
        plot (t, psolRedu012(n+1:2*n, :), LineWidth=2);
        hold on
        yline(0, '--');
        hold off
        xlabel('t');
        ylabel('Probability/\epsilon');
        fontsize(f2, scale=1.5);
        exportgraphics(f2, '8stateoeps.pdf');
        
        f3 = figure;
        plot (t, psolRedu012(2*n+1:3*n, :), LineWidth=2);
        hold on
        yline(0,  '--');
        xlabel('t');
        ylabel('Probability/\epsilon^2');
        fontsize(f3, scale=1.5);
        exportgraphics(f3, '8stateoeps2.pdf');
    end

    if plotVerbose    
        figure
        plot(t,Q*psolFull)
        set(gca,'xscale','log')
        title('Full Model - Projected Solution')
    end

    if plotVerbose
        figure
        plot(t,psolRedu0)
        set(gca,'xscale','log')
        title('O(1) Approximation')
    end

    if plotVerbose
        figure
        plot(t,psolRedu1)
        set(gca,'xscale','log')
        title('O(\epsilon) Approximation')

        figure
        plot(t,psolRedu012(n+1:2*n,:))
        set(gca,'xscale','log')
        title('O(\epsilon) Correction')
    end

    if plotVerbose
        figure
        plot(t,psolRedu2)
        set(gca,'xscale','log')
        title('O(\epsilon^2) Approximation')

        figure
        plot(t,psolRedu012(n+1:2*n,:))
        set(gca,'xscale','log')
        title('O(\epsilon) Correction')

        figure
        plot(t,psolRedu012(2*n+1:3*n,:))
        set(gca,'xscale','log')
        title('O(\epsilon^2) Correction')
    end

    if fairplot
        f = figure
        fontsize(f, 20, "points")
        C = colororder;
        colororder(C(1:4, :));
        a = plot(t, Q*psolFull, 'LineWidth', 2);
        hold on
        b = plot(t, psolRedu0, '--', 'LineWidth', 2);
        c = plot(t, psolRedu1, '-.', 'LineWidth', 2);
        d = plot(t, psolRedu2, ':', 'LineWidth', 2);
        hold off
        xlabel('t');
        ylabel('Probability');
        hleglines = [a(3), b(3), c(3), d(3)];
        % legend(a, 'Superbasin 1', 'Superbasin 2', 'Superbasin 3', 'Superbasin 4');
        % legend((a(1), b(1), c(1), d(1)),'Full Solution', ...
        %     'O(1) approximation', 'O(\epsilon) approximation', ...
        %     'O(\epsilon^2) approximation');
        % legend(a(1), b(1), 'Full Solution', 'O(1) approximation')
        legend(hleglines, 'Full solution', 'O(1) approximation', ...
            'O(\epsilon) approximation', 'O(\epsilon^2) approximation');
        ylim([0, 1]);
        title('\epsilon=0.005')
        fontsize(f, scale=1.5);
        % set(gca,'xscale','log')
        exportgraphics(f, 'basins0005.pdf')
    end

end

% close all

%% Plot figure over range of epsilon values

fff = figure
epsl_range
normErr0
normErr1
normErr2

save('8state.mat', 'epsl_range', 'normErr0', 'normErr1', 'normErr2')
% save('8stateuneq.mat', 'epsl_range', 'normErr0', 'normErr1', 'normErr2')
% save('8statenoneq.mat', 'epsl_range', 'normErr0', 'normErr1', 'normErr2')
% save('8stateeqo1.mat', 'epsl_range', 'normErr0', 'normErr1', 'normErr2')
% save('3state.mat', 'epsl_range', 'normErr0', 'normErr1', 'normErr2')

loglog(epsl_range,normErr0,'ob','MarkerFaceColor','b')
hold on
loglog(epsl_range,normErr1,'sr','MarkerFaceColor','r')
loglog(epsl_range,normErr2,'dk','MarkerFaceColor','k')
plot(epsl_range,0.21*epsl_range.^1,'-b')
plot(epsl_range,0.17*epsl_range.^2,'-r')
plot(epsl_range,0.28*epsl_range.^3,'-k')
% plot(epsl_range,0.0005*epsl_range.^1,'-b')
% plot(epsl_range,0.0007*epsl_range.^2,'-r')
% plot(epsl_range,0.0024*epsl_range.^3,'-k')
hold off
if length(epsl_range) > 1
    set(gca,'xlim',[min(epsl_range) max(epsl_range)])
end
hleg = legend('O(1) Approximation Error',...
    'O(\epsilon) Approximation Error',...
    'O(\epsilon^2) Approximation Error',...
    'C\times\epsilon',...
    'C\times\epsilon^2',...
    'C\times\epsilon^3');
set(hleg,'Location','SouthEast')
xlabel('\epsilon')
ylabel('Error')

print -dpng -r150 ApproxErrorVsEpsilon
