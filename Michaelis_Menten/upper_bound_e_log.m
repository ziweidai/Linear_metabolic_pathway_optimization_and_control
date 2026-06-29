function [c,ceq,gc,gceq] = upper_bound_e_log(loge)
e = 10.^loge;
c = 0;
ceq = sum(e)-1;
gc = 0;
gceq = log(10)*e;
end
