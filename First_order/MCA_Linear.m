function [S,f,fcc,dg]=MCA_Linear(a,e,K,Sin,Sout)
k = a.*e;
[S,f]=SteadyStateConc(k,K,Sin,Sout);
nrxn=length(k);
fcc = zeros(nrxn,1);
for i=1:nrxn
    fcc(i)=(S(i)*K(i)-S(i+1))/prod(K(1:i));
end
fcc=fcc/sum(fcc);
Q=S(2:end)./S(1:end-1);
dg=log(Q./K);
end
