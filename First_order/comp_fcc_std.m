function gstd = comp_fcc_std(a,e,K,Sin,Sout)
% Compute standard deviation of FCCs directly from pathway kinetic
% parameters
[~,~,fcc,~]=MCA_Linear(a,e,K,Sin,Sout);
gstd = std(fcc);
end