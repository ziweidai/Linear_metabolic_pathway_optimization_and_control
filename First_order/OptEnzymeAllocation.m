function e = OptEnzymeAllocation(a, K)
% K_i = (a_i*e_i^2)/(a_{i+1}*e_{i+1}^2)
% This means that e_i/e_{i+1} = sqrt(K_i*a_{i+1}/a_i), let e_ratio
% represent the right term
nrxn = length(a);
e_ratio = sqrt(K(1:nrxn-1).*a(2:nrxn)./a(1:nrxn-1));
e = ones(nrxn,1);
for i = nrxn-1:-1:1
    e(i) = e(i+1)*e_ratio(i);
end
e = e./sum(e);
