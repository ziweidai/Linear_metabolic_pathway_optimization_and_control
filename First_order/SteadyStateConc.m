function [S,f]=SteadyStateConc(k,K,Sin,Sout)
%Calculate steady state concentrations of metabolites
%k - rate coefficients (enzyme activities) of each step
%K - equilibrium constant of each step
%Sin - input substrate (nutrient) concentration
%Sout - output substrate (end product) concentration
nrxn=length(k);
nmet=nrxn+1;
b=zeros(nmet,1);
b(1)=Sin;
b(end)=Sout;

% Perform pre-conditioning of the linear equations to avoid numerical instability
logP = zeros(nmet,1);            
logP(2:end-1) = cumsum(log(K(2:end)));  
P = exp(logP - mean(logP));   

%{
M = zeros(nmet,nmet);
M(1,1)=1;
M(end,end)=1;
for i=1:nrxn-1
    M(i+1,i)   = k(i);
    M(i+1,i+2) = k(i+1)/K(i+1);
    M(i+1,i+1) = -(k(i)/K(i) + k(i+1));
end
A = M * diag(P);
%}

dl = zeros(nmet,1); d = zeros(nmet,1); du = zeros(nmet,1);
d(1)=1; d(end)=1;
for i=1:nrxn-1
    dl(i)=k(i);
    du(i+2)=k(i+1)/K(i+1);
    d(i+1)=-(k(i)/K(i) + k(i+1));
end
M = spdiags([dl d du],[-1 0 1],nmet,nmet);
A = M * spdiags(P,0,nmet,nmet);


u = A \ b;
S = P .* u;

f=k(1)*S(1)-k(1)*S(2)/K(1);
end