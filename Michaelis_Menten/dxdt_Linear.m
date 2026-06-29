function [dxdt,J]=dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq)
%Calculate dxdt for linear pathway with Michaelis-Menten kinetics
n = length(x);
S = [Sin;x];
P = [x;Sout];
v=MM(S,P,kcat,Ks,Kp,Keq);
dxdt = v(1:n)-v(2:n+1); %net changes of metabolites
J=Jacobian_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq);
end
