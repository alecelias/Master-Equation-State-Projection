clear all
% close all
clc

global n W S Q QT FNorm invFNorm indxpiv epsl

%% Set up the problem
% Define the number of states and the propensities for the allowed transitions
% kfij refers to transition from state i to j ( i => j )

% RunParSet1
% RunParSet2
% RunParSet3
% disp(['Number of states: ' num2str(m)]);

% epsl_range = [3.2 1.6 0.8 0.4 0.2 0.1 0.05 0.025 0.0125 0.00625 0.003125];
% epsl_range = logspace(-2.1,-4.1,30);
epsl_range = [0.01];
normErr0 = [];
normErr1 = [];
normErr2 = [];
maxNormalizErrFullSys = [];
maxNormalizErrReduSys0 = [];
maxNormalizErrReduSys1 = [];
maxNormalizErrReduSys2 = [];
avesss = [];

plotVerbose = false;
% plotVerbose = true;

% plotave = false;
plotave = true;

snapplot = false;
% snapplot = true;

fairplot = false;
% fairplot = true;

% Set up the transition matrix for 1D lattice

% L = 9;
% m = 2^L;
% [F, S] = CreateTmats(L, 1, 8.5, 4.0, 10.);

% Set up the transition matrix for 2D lattice

L = 3;
NL = L^2;
m = 2^NL;
[F, S] = Create2Dmats(L, 1, 8.5, 4.0, 10.);

%% Create necessary matrices

Nsp = round(rref(transpose(null(transpose(F))))); % left nullspace of F
n = rank(Nsp); % number of superbasins
Q = Nsp;
% QT = QT(:,n:-1:1); % rearrange columns (for convenient relabelling of superbasins)
% QT = QT(:,[3 1 4 2]); % some arbitrary permutation of columns
QT = transpose(Q);
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
% ODESolOptions = odeset('RelTol',1e-7,'AbsTol',1e-10); % Options for the ODE solver
ODESolOptions = odeset('RelTol',1e-12,'AbsTol',1e-14); % Options for the ODE solver
tspan = [0 1];
t = linspace(tspan(1),tspan(2),1001); % times to evaluate the solution
Dt = t(2)-t(1);
%% Create initial conditions
p = [0.15, 0.05, 0.08, 0.12, 0.16, 0.04, 0.13, 0.07, 0.03, 0.17]';
rhstmp = QT*p;
% p = rand(m, 1);
% p = p / sum(p);
% rhstmp = QT*Q*p;
rhstmp(~indxpiv) = 0;
pini0 = invFNorm*rhstmp; % use quasi equilibrated (within superbasins) probabilities

ptilde012ini = [Q*pini0; zeros(2*n,1)];

%% Solve the approximate system
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

    if plotVerbose
        figure
        plot(t,psolFull)
        set(gca,'xscale','log')
        title('Full Model - Full Solution')

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
        plot(t,psolRedu01(n+1:2*n,:))
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

    if plotave
        rFull = Q * psolFull;
        nads = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        aveFull = nads * rFull;
        avesss = [avesss; aveFull];
        fave = figure;
        ave0 = nads * psolRedu0;
        ave1 = nads * psolRedu1;
        ave2 = nads * psolRedu2;
        plot(t, aveFull, LineWidth=3);
        hold on
        plot(t, ave0, '--', LineWidth=2);
        plot(t, ave1, '-.', LineWidth=1.5);
        plot(t, ave2, ':', LineWidth=1);
        xlabel('t');
        ylabel('<N_{ads}>');
        legend('Full Solution', 'O(1) approximation', 'O(\epsilon) approximation', 'O(\epsilon^2) approximation');
        % ylabel('\left< N_{ads} \right>', 'Interpreter','latex');
        % hold on
        % plot(t, ave0, '--');
        % plot(t, ave1, '-.');
        % plot(t, ave2, ':');
        % hold off
        ylim([2.5, 4.5]);
        hold off
        set(gca,'xscale','log');
        title('\epsilon = 0.01');
        fontsize(fave, scale=1.5);
        exportgraphics(fave, 'latticeaverage.pdf')

    end

    if snapplot
        rFull = Q * psolFull;
        f1 = figure;
        nads = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        a = bar(nads, rFull(:, 1));
        hold on
        c = scatter(nads - 0.2222, psolRedu0(:, 1), 'LineWidth', 2);
        d = scatter(nads, psolRedu1(:, 1), 'x', 'LineWidth', 2);
        e = scatter(nads + 0.2222, psolRedu2(:, 1), 'd', 'LineWidth', 2);
        b = bar(nads, [psolRedu0(:, 1), psolRedu1(:, 1), psolRedu2(:, 1)], 0);
        xlabel('N_{ads}');
        ylabel('Probability');
        hold off
        fontsize(f1, scale=1.5);
        ylim([0, 0.5]);
        title('\epsilon = 0.01, t = 0');
        legend('Full solution', 'O(1) approximation', 'O(\epsilon) approximation', 'O(\epsilon^2) approximation');
        exportgraphics(f1, 'latticeinitial.pdf');

        f2 = figure;
        bar(nads, rFull(:, 751));
        hold on
        c = scatter(nads - 0.2222, psolRedu0(:, 751), 'LineWidth', 2);
        d = scatter(nads, psolRedu1(:, 751), 'x', 'LineWidth', 2);
        e = scatter(nads + 0.2222, psolRedu2(:, 751), 'd', 'LineWidth', 2);
        b = bar(nads, [psolRedu0(:, 751), psolRedu1(:, 751), psolRedu2(:, 751)], 0);
        xlabel('N_{ads}');
        ylabel('Probability');
        hold off
        legend('Full solution', 'O(1) approximation', 'O(\epsilon) approximation', 'O(\epsilon^2) approximation');
        fontsize(f2, scale=1.5);
        ylim([0, 0.5]);
        title('\epsilon = 0.01, t = 0.75');
        exportgraphics(f2, 'latticet750.pdf');

        f3 = figure;
        bar(nads, rFull(:, 31));
        hold on
        c = scatter(nads - 0.2222, psolRedu0(:, 31), 'LineWidth', 2);
        d = scatter(nads, psolRedu1(:, 31), 'x', 'LineWidth', 2);
        e = scatter(nads + 0.2222, psolRedu2(:, 31), 'd', 'LineWidth', 2);
        b = bar(nads, [psolRedu0(:, 31), psolRedu1(:, 31), psolRedu2(:, 31)], 0);
        xlabel('N_{ads}');
        ylabel('Probability');
        hold off
        fontsize(f3, scale=1.5);
        ylim([0, 0.5]);
        title('\epsilon = 0.01, t = 0.03');
        legend('Full solution', 'O(1) approximation', 'O(\epsilon) approximation', 'O(\epsilon^2) approximation');
        exportgraphics(f3, 'latticet030.pdf');
    end

    if fairplot
        f = figure
        fontsize(f, 20, "points")
        C = colororder;
        % colororder(C(1:4, :));
        a = plot(t, Q*psolFull, 'LineWidth', 2);
        % hold on
        % b = plot(t, psolRedu0, '--', 'LineWidth', 2);
        % c = plot(t, psolRedu1, '-.', 'LineWidth', 2);
        % d = plot(t, psolRedu2, ':', 'LineWidth', 2);
        % hold off
        xlabel('t');
        ylabel('Probability');
        set(gca,'xscale','log')

        % hleglines = [a(3), b(3), c(3), d(3)];
        % legend(a, 'Superbasin 1', 'Superbasin 2', 'Superbasin 3', 'Superbasin 4');
        % legend((a(1), b(1), c(1), d(1)),'Full Solution', ...
        %     'O(1) approximation', 'O(\epsilon) approximation', ...
        %     'O(\epsilon^2) approximation');
        % legend(a(1), b(1), 'Full Solution', 'O(1) approximation')
        % legend(hleglines, 'Full solution', 'O(1) approximation', ...
        %     'O(\epsilon) approximation', 'O(\epsilon^2) approximation');
        % ylim([0, 1]);
    end
end

% close all
figure
epsl_range
normErr0
normErr1
normErr2

% save('1Dlat.mat', 'epsl_range', 'normErr0', 'normErr1', 'normErr2')
% save('2Dlat.mat', 'epsl_range', 'normErr0', 'normErr1', 'normErr2')

loglog(epsl_range,normErr0,'ob','MarkerFaceColor','b')
hold on
loglog(epsl_range,normErr1,'sr','MarkerFaceColor','r')
loglog(epsl_range,normErr2,'dk','MarkerFaceColor','k')
% plot(epsl_range,0.21*epsl_range.^1,'-b')
% plot(epsl_range,0.17*epsl_range.^2,'-r')
% plot(epsl_range,0.28*epsl_range.^3,'-k')
% plot(epsl_range,0.0005*epsl_range.^1,'-b')
% plot(epsl_range,0.0007*epsl_range.^2,'-r')
% plot(epsl_range,0.0024*epsl_range.^3,'-k')
plot(epsl_range,0.20*epsl_range.^1,'-b')
plot(epsl_range,1.14*epsl_range.^2,'-r')
plot(epsl_range,6.27*epsl_range.^3,'-k')
% plot(epsl_range,0.56*epsl_range.^1,'-b')
% plot(epsl_range,18.36*epsl_range.^2,'-r')
% plot(epsl_range,707.1*epsl_range.^3,'-k')
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

% if plotave
%     plot(t, avesss);
% end

% print -dpng -r150 ApproxErrorVsEpsilon
