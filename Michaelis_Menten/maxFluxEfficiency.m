function [J,e,Jflag] = maxFluxEfficiency(Sin,Sout,kcat,Ks,Kp,Keq)
nrxn = length(kcat);
x0 = ones(nrxn,1)/nrxn;
A = ones(1,nrxn);
b = 1;
lb = zeros(nrxn,1);
ub = ones(nrxn,1);

%options = optimoptions(@fmincon,'MaxIterations',5000, ...
%    'OptimalityTolerance',1e-10);
options = optimoptions('fmincon',...
    'Display','none',...       
    'Algorithm','interior-point',...
    'OptimalityTolerance',1e-6,...
    'ConstraintTolerance',1e-6,...
    'StepTolerance',1e-10,...
    'MaxIterations',5000,...
    'MaxFunctionEvaluations',1e5);

try
    [e,~,exitflag] = fmincon(@(x)-SS_Linear(Sin,Sout,kcat.*x,Ks,Kp,Keq), ...
        x0,[],[],A,b,lb,ub,[],options);
    [J,~,Jflag] = SS_Linear(Sin,Sout,kcat.*e,Ks,Kp,Keq);
    if exitflag <= 0 % filter out solutions with failed optimization
        warning('Optimization failed with exitflag = %d', exitflag)
        e = x0;
        J = -1;
    end
catch ME
    warning('An error occurred within fmincon: %s', ME.message);
    e = x0;
    J = -1;
    Jflag = -3;
end
end