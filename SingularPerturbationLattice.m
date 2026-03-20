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
epsl_range = logspace(0.5,-3,30);
% epsl_range = [0.01];
normErr0 = [];
normErr1 = [];
normErr2 = [];
maxNormalizErrFullSys = [];
maxNormalizErrReduSys0 = [];
maxNormalizErrReduSys1 = [];
maxNormalizErrReduSys2 = [];

plotVerbose = false;
% plotVerbose = true;

for epsl = epsl_range

    [F, S] = CreateTmats(4, 1E+9, 8.5E+5, 4.0E+5, 1E+6)
    m = 2^4
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

    p = zeros(m, 1);
    p(2) = 1;
    % Demonstrate the solution of the O(1/eps) equations
    rhstmp = QTQ*p; % right-hand side must have the values of ptilde
    % in the rows corresponding to the normalisation equations
    rhstmp(~indxpiv) = 0; % the other rows should have zeros (fluxes between superbasin's states balance out)
    pQEq = invFNorm*rhstmp; % this gives the probabilities quasi-equilibrated within each superbasin
    sum(pQEq);
    % Validate the projector: the following expression should evaluate to the
    % zero matrix
    Q*F;

    % Check the calculation of conditional probabilities
    p = zeros(m, 1);
    p(2) = 1;
    % p = [0.3 0.5 0.2].';

    % Check the validity of the distribution matrix K
    K = diag(p)*inv(diag(QT*Q*p))*QT;
    ptilde = Q*p;
    p - K*ptilde; % this should give the zero vector (up to numerical accuracy)


    %% Solve the full and reduced models
    % ODESolOptions = odeset('RelTol',1e-7,'AbsTol',1e-10); % Options for the ODE solver
    ODESolOptions = odeset('RelTol',1e-10,'AbsTol',1e-12); % Options for the ODE solver
    tspan = [0 10];
    t = linspace(tspan(1),tspan(2),1001); % times to evaluate the solution
    Dt = t(2)-t(1);

    %% Create initial conditions (this is still experimental;
    %% it is unclear whether it provides any benefits)
    rhstmp = QT*Q*p;
    rhstmp(~indxpiv) = 0;
    pini0 = invFNorm*rhstmp; % use quasi equilibrated (within superbasins) probabilities

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

    % Solve the full system
    solFull = ode15s(@FullSystemRHS,tspan,pini,ODESolOptions);
    psolFull = deval(solFull,t);
    maxNormalizErrFullSys = [maxNormalizErrFullSys max(sum(psolFull,1)-1)];

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

    % Solve the O(1) reduced system
    ptilde0ini = Q*pini;
    solRedu0 = ode15s(@ReduOe0SystemRHS,tspan,ptilde0ini,ODESolOptions);
    psolRedu0 = deval(solRedu0,t);
    maxNormalizErrReduSys0 = [maxNormalizErrReduSys0 max(sum(psolRedu0,1)-1)];

    if plotVerbose
        figure
        plot(t,psolRedu0)
        set(gca,'xscale','log')
        title('O(1) Approximation')
    end

    % Evaluate the error of the O(1) system
    Err0 = sqrt(sum((Q*psolFull-psolRedu0).^2*Dt,2));
    normErr0 = [normErr0 norm(Err0)];

    % Solve the O(eps) reduced system
    ptilde01ini = [Q*pini; zeros(n,1)];
    solRedu01 = ode15s(@ReduOe1SystemRHS,tspan,ptilde01ini,ODESolOptions);
    psolRedu01 = deval(solRedu01,t);
    psolRedu1 = psolRedu01(1:n,:) + epsl*psolRedu01(n+1:2*n,:);
    maxNormalizErrReduSys1 = [maxNormalizErrReduSys1 max(sum(psolRedu1,1)-1)];

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

    % Evaluate the error of the O(eps) system
    Err1 = sqrt(sum((Q*psolFull-psolRedu1).^2*Dt,2));
    normErr1 = [normErr1 norm(Err1)];

    % Solve the O(eps^2) reduced system
    % ptilde012ini = [Q*pini; zeros(n,1); [-0.003226733393321; 0.003226733393321]];
    ptilde012ini = [Q*pini; zeros(2*n,1)];
    solRedu012 = ode15s(@ReduOe2SystemRHS,tspan,ptilde012ini,ODESolOptions);
    psolRedu012 = deval(solRedu012,t);
    psolRedu2 = psolRedu012(1:n,:) + epsl*psolRedu012(n+1:2*n,:) ...
        + epsl^2*psolRedu012(2*n+1:3*n,:);
    maxNormalizErrReduSys2 = [maxNormalizErrReduSys2 max(sum(psolRedu2,1)-1)];

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

    % Evaluate the error of the O(eps^2) system
    Err2 = sqrt(sum((Q*psolFull-psolRedu2).^2*Dt,2));
    normErr2 = [normErr2 norm(Err2)];

end

% close all
figure
epsl_range
normErr0
normErr1
normErr2

save('1Dlat.mat', 'epsl_range', 'normErr0', 'normErr1', 'normErr2')

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
