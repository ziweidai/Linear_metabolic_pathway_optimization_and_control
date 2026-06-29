met_conc_table = readtable("Metabolite_concentration_rabinowitz_2016.csv");
kcat_Km_table = readtable("Human_kcat_Km.xlsx");
deltaG_table = readtable("Recon3D_standard_dGr_dGbyG.csv");

met_conc_table = unique(met_conc_table,"rows");
kcat_Km_table = unique(kcat_Km_table,"rows");
deltaG_table = unique(deltaG_table,"rows");
met_conc_table.Organism = string(met_conc_table.Organism);
met_conc_table = met_conc_table(met_conc_table.Organism == "Homo sapiens", :);

met_conc_data = met_conc_table.Concentration_M_;
kcat_data = kcat_Km_table.Kcat_predicted_;
Km_data = kcat_Km_table.KM_mM__predicted_/1000; %Convert unit to mol/L
deltaG_data = deltaG_table.standard_dGr_prime_kJ_mol_*1000;
R = 8.314;
T = 298.15;
K = exp(-deltaG_data/R/T);
K = K(~isnan(K));
K = K(K>1e-10 & K<1e10);

writematrix(K,"K.csv");
writematrix(kcat_data,"kcat.csv");
writematrix(Km_data,"Km.csv");
writematrix(met_conc_data,"S.csv");