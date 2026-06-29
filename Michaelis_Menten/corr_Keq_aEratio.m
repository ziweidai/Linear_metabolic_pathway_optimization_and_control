function [c,c_adj,r,r_adj,sat,phi,ae2r] = corr_Keq_aEratio(Keq,kcat,Ks,Kp,Sin,Sout,e)
%Compute the correlation between Keq and the ratio of aE^2
%a is computed as kcat/(1+S/Ks+P/Kp)
%sat is the saturation term, computed by the function Saturation[]
n = length(kcat) - 1;
y0=log10(((Sin*Sout)^0.5)*ones(n,1));
options = optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-10);
%y=fsolve(@(x) dxdt_Linear(x,Sin,Sout,kcat,Ks,Kp,Keq),y0,options);
kcat_appr = kcat.*e;
[y,~] = lsqnonlin(@(x) dxdt_Linear(10.^x,Sin,Sout,kcat_appr,Ks,Kp,Keq), ...
    y0,[],[],options);
y = 10.^y;
S_ext=[Sin;y(:);Sout];
S = S_ext(1:end-1);
P = S_ext(2:end);
a = kcat./Ks./(1+S./Ks+P./Kp);
ae2 = a.*e.*e;
c = corr(log10(Keq(1:end-1)),log10(ae2(1:end-1)./ae2(2:end)));
r = rmse(log10(Keq(1:end-1)),log10(ae2(1:end-1)./ae2(2:end)));
sat = 1-1./(1 + S./Ks + P./Kp);
phiS = (Kp.*P + Ks.*Keq.*(P+Kp))./(Ks.*Keq.*P + Kp.*Keq.*(Ks+S));
phiP = (Ks.*Keq.*S + Kp.*(Ks+S))./(Ks.*P + Kp.*(Ks+S));
phi = phiS(2:n+1)./phiP(1:n);
c_adj = corr(log10(Keq(1:n).*phi),log10(ae2(1:end-1)./ae2(2:end)));
r_adj = rmse(log10(Keq(1:end-1).*phi),log10(ae2(1:end-1)./ae2(2:end)));
ae2r = ae2(1:end-1)./ae2(2:end);
%{
[c c_adj]
v1 = MM(S(1),P(1),kcat(1)*e(1),Ks(1),Kp(1),Keq(1));
v2 = MM(S(2),P(2),kcat(2)*e(2),Ks(2),Kp(2),Keq(2));
if v1>0 && abs(v1-v2)<1e-4
    figure;
    subplot(1,2,1);scatter(log(Keq(1:n)),log(ae2(1:end-1)./ae2(2:end)));
    subplot(1,2,2);scatter(log(Keq(1:n).*phi),log(ae2(1:end-1)./ae2(2:end)));
end
%}
end
