function [phiS,phiP]=calculateRxnPhi(S,P,Ks,Kp,Keq)
phiS = (Kp*P + Ks*Keq*(P+Kp))/(Ks*Keq*P + Kp*Keq*(Ks+S));
phiP = (Ks*Keq*S + Kp*(Ks+S))/(Ks*P + Kp*(Ks+S));
end