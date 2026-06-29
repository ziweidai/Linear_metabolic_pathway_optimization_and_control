%% Analyze distribution of w_i for glycolysis
% Read thermodynamic and kinetic parameters for glycolysis from data files
species = "ecoli"; %Specify the species (human, yeast, or ecoli)
pathway = "TCA"; %Specify the pathway (glycolysis or TCA, or glycolysis_with_GLUT)

filename_K = sprintf("data/%s_K_%s.csv",pathway,species);
filename_kcat_Km = sprintf("data/%s_catalytic_efficiency_%s.csv",pathway,species);
filename_met_conc = sprintf("data/%s_met_concentration.csv",pathway);
filename_CAQ = sprintf("data/%s_%s_mean_CAQ.csv",species,pathway);

K_table = readtable(filename_K,"ReadVariableNames",true);
met_conc_table = readtable(filename_met_conc, ...
    "ReadVariableNames",true,"ReadRowNames",true);
met_conc_exp = met_conc_table{:,species};
CAQ_exp = readtable(filename_CAQ);

met_names = met_conc_table.metname;
rxn_names = K_table.rxnname;
rxn_pair_names = rxn_names(1:end-1);
for i = 1:length(rxn_names)-1
    rxn_pair_names{i} = sprintf("%s-%s",rxn_names{i},rxn_names{i+1});
end

%Correct the equilibrium constant for ALDO
K_pathway = K_table.K;

if strcmp(pathway,"glycolysis")
    K_pathway(4) = K_pathway(4)/(met_conc_table{"g3p",species}/1000);
elseif strcmp(pathway,"glycolysis_with_GLUT")
    K_pathway(5) = K_pathway(5)/(met_conc_table{"g3p",species}/1000);
end

kcat_Km_table = readtable(filename_kcat_Km, ...
    "ReadVariableNames",true);

n_rxn_pathway = size(K_table,1);
a_pathway = ones(n_rxn_pathway,1);
for i = 1:n_rxn_pathway
    rxn_name = K_table.Reaction(i);
    rxn_idx = find(ismember(kcat_Km_table.BiGGReactionID,rxn_name));    
    a0 = kcat_Km_table.enzymeEffeciency(rxn_idx);
    enzyme_mass = kcat_Km_table.enzymeMass(rxn_idx);
    a = a0./enzyme_mass;
    if length(rxn_idx)>1
        a_pathway(i) = geomean(a);
    else
        a_pathway(i) = a;
    end
end

% Set model parameters by the pathway parameters
nrxn = n_rxn_pathway;
nmet = nrxn-1;
K = K_pathway;
a = a_pathway;

if strcmp(pathway,"glycolysis") 
    Sin = met_conc_table{"glc",species};
    Sout = met_conc_table{"pyr",species}/2;
elseif strcmp(pathway,"glycolysis_with_GLUT")
    Sin = met_conc_table{"glc_D_EX",species};
    Sout = met_conc_table{"pyr",species}/2;
elseif strcmp(pathway,"TCA")
    Sin = met_conc_table{"cit",species};
    Sout = met_conc_table{"oaa",species};
end

if strcmp(species,"human")
    species_name = "H. sapiens";
elseif strcmp(species,"yeast")
    species_name = "S. cerevisiae";
else
    species_name = "E. coli";
end

if strcmp(pathway,"TCA")
    pathway_name = "TCA cycle";
else
    pathway_name = "Glycolysis";
end
condition_label = sprintf("%s, %s", pathway_name, species_name);

%% Simulate random solutions, compare flux/control efficiency
nsamples = 2000;
eSample = 10.^(4*lhsdesign(nsamples,nrxn)');
%normalize eSample so that sum of enzyme abundances equals 1
eSample = eSample./sum(eSample);

SMat = zeros(nsamples,nmet+2);
fVec = zeros(nsamples,1);
fccMat = zeros(nsamples,nrxn);

for i = 1:nsamples
    [S,f,fcc,dg] = MCA_Linear(a, eSample(:,i), K, Sin, Sout);
    SMat(i,:) = S';
    fVec(i) = f;
    fccMat(i,:)=fcc';
end
SMat = SMat';
SMat = SMat(:,fVec>0);
SMat = adjust_met_conc(SMat,pathway);
e_max_flux = OptEnzymeAllocation(a, K);
[S_max_flux,f_max,fcc_max_flux,~] = MCA_Linear(a, e_max_flux, K, Sin,Sout);
fcc_std_max_flux = std(fcc_max_flux);
S_max_flux = adjust_met_conc(S_max_flux,pathway);

%% Simulate the Pareto front under varying flux
log10_flux_scale = -0.03:-0.03:-3;
coef_red_flux = 10.^(log10_flux_scale); % Varying flux from 0.001f_max to f_max
n_points = length(coef_red_flux);
Pareto_fcc_std_values = zeros(n_points,1);
e_max_std_varying_f = zeros(nrxn,n_points);
S_std_varying_f = zeros(nrxn+1, n_points);
FCC_varying_f = zeros(nrxn,n_points);
e0 = e_max_flux;
Aeq = ones(1,nrxn);
beq = 1;
lb = 1e-4*ones(nrxn,1);
ub = ones(nrxn,1);
options = optimoptions('fmincon',...
    'Display','none',...       
    'Algorithm','interior-point',...
    'OptimalityTolerance',1e-6,...
    'ConstraintTolerance',1e-6,...
    'StepTolerance',1e-10,...
    'MaxIterations',5000,...
    'MaxFunctionEvaluations',1e5);
exitflag_std = zeros(n_points,1);
firstorder_std = zeros(n_points,1);
constrviol_std = zeros(n_points,1);
for i = 1:n_points
    f_min = coef_red_flux(i)*f_max;
    
    if i>1
        e0 = e_max_std_varying_f(:,i-1);
    end
    [e_max_std_varying_f(:,i),Pareto_fcc_std_values(i),...
        exitflag_std(i),output] = ...
        fmincon(@(e)-comp_fcc_std(a,e,K,Sin,Sout),e0,[],[], ...
        Aeq,beq,lb,ub,@(e)flux_bound(a,e,K,Sin,Sout,f_min),options);
    [S_std_varying_f(:,i),~] = SteadyStateConc(a.*e_max_std_varying_f(:,i),K,Sin,Sout);
    [S,f,fcc,dg] = MCA_Linear(a, e_max_std_varying_f(:,i), K, Sin, Sout);
    FCC_varying_f(:,i) = fcc;
        
    firstorder_std(i) = output.firstorderopt;
    constrviol_std(i) = output.constrviolation;
end
Pareto_fcc_std_values = -Pareto_fcc_std_values;
idx_std_good = find(exitflag_std>0 & firstorder_std < 1e-4 & constrviol_std < 1e-4);
S_std_varying_f = adjust_met_conc(S_std_varying_f,pathway);
fcc_std = std(fccMat,[],2);


%% Compare the real and model-predicted [S] (the one for Pareto solution closest to real CAQ)
caq_max_std_varying_f = a(1:end-1).*(e_max_std_varying_f(1:end-1,:).^2)...
    ./(a(2:end).*(e_max_std_varying_f(2:end,:).^2));
log10CAQ_exp = log10(table2array(CAQ_exp));
log10CAQ_Pareto = log10(caq_max_std_varying_f);
dev_CAQ = sqrt(sum((log10CAQ_exp(:,1)-log10CAQ_Pareto).^2));
[~,idx_min_dev_caq] = min(dev_CAQ);
figure;
x = S_std_varying_f(:,idx_min_dev_caq);
y = met_conc_exp;

rmse_S = rmse(log10(x),log10(y));
pearsonR_S = corr(log10(x),log10(y));

blues = brewermap(2,'Blues');
blues = blues(end:-1:1,:);
log10range = [fix(log10(min([x;y])))-1 fix(log10(max([x;y])))+1];
range = 10.^log10range;
hold on;
h1 = plot(range,range,':','Color',blues(1,:));
h2 = plot(range,range*10,':','Color',blues(2,:));
plot(range,range/10,':','Color',blues(2,:));

scatter(x,y,64,'Marker','diamond','MarkerEdgeColor',[233 113 50]/256, ...
    'MarkerFaceColor',[233 113 50]/256);
annotation('textbox',[0.6 0.3 0.6 0.1], ...
    'String',sprintf("RMSE=%.1f\nPearson's R=%.2f",rmse_S,pearsonR_S), ...
    'EdgeColor','none', 'FontSize',20);

automatic_label(x, y, met_names, 1, [0.2 0.02]);


xlabel("[S], CAQ-matched Pareto solution");
ylabel("[S], experimental");
title(condition_label);
set(gca,'xscale','log');
set(gca,'yscale','log');
legend([h1, h2], {'Equal', '10-fold difference'}, ...
    'Location', 'best','FontSize',20);
legend boxoff;
xlim(range);
ylim(range);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_S_Pareto_matched.fig",condition_label);
pdf_name = sprintf("%s_S_Pareto_matched.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

data_output = table(S_std_varying_f(:,idx_min_dev_caq),S_max_flux,met_conc_exp);
data_output.Properties.RowNames = met_names;
data_output.Properties.VariableNames = ["CAQ-matched","Maximal_flux","Experimental"];
writetable(data_output,sprintf("%s_metconc_comparison.csv",condition_label),...
    'WriteRowNames',true);

%% Compare the real and model-predicted [S] under maximal flux
figure;
x = S_max_flux;
y = met_conc_exp;

rmse_S = rmse(log10(x),log10(y));
pearsonR_S = corr(log10(x),log10(y));

blues = brewermap(2,'Blues');
blues = blues(end:-1:1,:);
log10range = [fix(log10(min([x;y])))-1 fix(log10(max([x;y])))+1];
range = 10.^log10range;
hold on;
h1 = plot(range,range,':','Color',blues(1,:));
h2 = plot(range,range*10,':','Color',blues(2,:));
plot(range,range/10,':','Color',blues(2,:));

scatter(x,y,64,'Marker','diamond','MarkerEdgeColor',[233 113 50]/256, ...
    'MarkerFaceColor',[233 113 50]/256);
annotation('textbox',[0.6 0.3 0.6 0.1], ...
    'String',sprintf("RMSE=%.1f\nPearson's R=%.2f",rmse_S,pearsonR_S), ...
    'EdgeColor','none', 'FontSize',20);

automatic_label(x, y, met_names, 1, [0.2 0.02]);

xlabel("[S], maximal flux");
ylabel("[S], experimental");
title(condition_label);
set(gca,'xscale','log');
set(gca,'yscale','log');
legend([h1, h2], {'Equal', '10-fold difference'}, ...
    'Location', 'best','FontSize',20);
legend boxoff;
xlim(range);
ylim(range);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_S_max_flux.fig",condition_label);
pdf_name = sprintf("%s_S_max_flux.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

%% Compare the real and model-predicted CAQ (closest to real)
figure;
x = caq_max_std_varying_f(:,idx_min_dev_caq);
y = table2array(CAQ_exp(:,1));
log10range = [fix(log10(min([x;y])))-1 fix(log10(max([x;y])))+1];
range = 10.^log10range;
hold on;
h1 = plot(range,range,':','Color',blues(1,:));
h2 = plot(range,range*10,':','Color',blues(2,:));
plot(range,range/10,':','Color',blues(2,:));
scatter(x,y,64,'Marker','diamond','MarkerEdgeColor',[233 113 50]/256, ...
    'MarkerFaceColor',[233 113 50]/256);
automatic_label(x, y, rxn_pair_names, 4, [1.5 0.02]);

rmse_caq = rmse(log10(x),log10(y));
pearsonR_caq = corr(log10(x),log10(y));
annotation('textbox',[0.6 0.3 0.6 0.1], ...
    'String',sprintf("RMSE=%.1f\nPearson's R=%.2f",rmse_caq,pearsonR_caq), ...
    'EdgeColor','none','FontSize',20);

legend([h1, h2], {'Equal', '10-fold difference'}, ...
    'Location', 'best','FontSize',20);
legend boxoff;
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim(range);
ylim(range);
xlabel("CAQ, CAQ-matched Pareto solution");
ylabel("CAQ, experimental");
title(condition_label);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);

fig_name = sprintf("%s_caq_Pareto_matched.fig",condition_label);
pdf_name = sprintf("%s_caq_Pareto_matched.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

data_output = table(x,y);
data_output.Properties.RowNames = string(rxn_pair_names);
data_output.Properties.VariableNames = ["CAQ-matched","Experimental"];
writetable(data_output,sprintf("%s_CAQ_comparison.csv",condition_label),...
    'WriteRowNames',true);

%% Plot the K-CAQ relationship on the Pareto front
flux_color = log10_flux_scale;
flux_color_norm = (flux_color - min(flux_color))/(max(flux_color) - min(flux_color));
flux_color_idx = round(flux_color_norm * 255) + 1;
cmap = brewermap(256,'GnBu');
%cmap = cmap(end:-1:1,:);
figure;
hold on;
data_range = [min([caq_max_std_varying_f(:);K]) max([caq_max_std_varying_f(:);K])];
range = 10.^[fix(log10(data_range(1)))-1 fix(log10(data_range(2)))+1];
h_leg0 = plot(range,range,':','Color',[0.2 0.2 0.2]);
for i = 1:length(idx_std_good)
    scatter(caq_max_std_varying_f(:,idx_std_good(i)),K(1:end-1),'Marker','+', ...
        'MarkerEdgeColor',cmap(flux_color_idx(idx_std_good(i)),:), 'LineWidth', 1.2);
end

colormap(cmap);
cb = colorbar;
cb.Label.String = '\eta/\eta_{max}';
cb.Ticks = [0 1/3 2/3 1];
cb.TickLabels = {'0.001','0.01','0.1','1'};
xlabel("CAQ");
ylabel("K");
%title("K-CAQ relationship on the Pareto front");
title(condition_label);
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim(range);
ylim(range);

h_leg2 = scatter(caq_max_std_varying_f(:,idx_std_good(idx_min_dev_caq)),...
    K(1:end-1),64,'Marker','diamond','MarkerEdgeColor',[233 113 50]/256, ...
    'MarkerFaceColor',[233 113 50]/256);

automatic_label(caq_max_std_varying_f(:,idx_std_good(idx_min_dev_caq)), ...
    K(1:end-1), rxn_pair_names, 4, [1.5 0.02]);

data_output = array2table([K(1:end-1) caq_max_std_varying_f(:,idx_std_good)]);
data_output.Properties.RowNames = string(rxn_pair_names);
CAQ_labels = compose("CAQ_%d",1:length(idx_std_good));
data_output.Properties.VariableNames = ["K",CAQ_labels];
writetable(data_output,sprintf("%s_Pareto_K_CAQ.csv",condition_label),...
    'WriteRowNames',true);

% Generate dummy handle for legend
h_leg1 = plot(nan, nan, '+k', 'MarkerSize', 8, 'LineWidth', 1.2);
legend([h_leg1, h_leg2, h_leg0], ...
    {'Pareto-front', 'CAQ-matched Pareto solution', 'K=CAQ line'}, ...
    'Location', 'best','FontSize',20);
legend boxoff;
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_Pareto_K_CAQ.fig",condition_label);
pdf_name = sprintf("%s_Pareto_K_CAQ.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

%% Plot the Pareto front
figure;
scatter(fVec(fVec>0),fcc_std(fVec>0),64,'Marker','.','MarkerEdgeColor',[0.7 0.7 0.7]);

data_output = table(fVec(fVec>0),fcc_std(fVec>0));
data_output.Properties.VariableNames = ["Flux_efficiency","Control_efficiency"];
writetable(data_output,sprintf("%s_flux_vs_control_random.csv",condition_label));

hold on;
plot([f_max f_max*coef_red_flux(idx_std_good)],...
    [fcc_std_max_flux Pareto_fcc_std_values(idx_std_good)'], ...
    'Color',[78 167 46]/256,'LineWidth',4);

data_output = table([f_max f_max*coef_red_flux(idx_std_good)]',...
    [fcc_std_max_flux;Pareto_fcc_std_values(idx_std_good)]);
data_output.Properties.VariableNames = ["Flux_efficiency","Control_efficiency"];
writetable(data_output,sprintf("%s_flux_vs_control_Pareto.csv",condition_label));

scatter(f_max,fcc_std_max_flux,96,'Marker','diamond', ...
    'MarkerEdgeColor',[0 112 192]/256, 'MarkerFaceColor',[0 112 192]/256);
scatter(f_max*coef_red_flux(idx_std_good(idx_min_dev_caq)),...
    Pareto_fcc_std_values(idx_std_good(idx_min_dev_caq)), 96,...
    'Marker','diamond','MarkerEdgeColor',[233 113 50]/256,'MarkerFaceColor',[233 113 50]/256);
set(gca,'xscale','log');
ylim([0.25 0.4]);%yticks([0.315 0.316 0.317]);
ylabel("Control efficiency (\kappa)");
xlabel("Flux efficiency (\eta)");
title(condition_label);
legend("Random","Pareto front","Maximal flux","CAQ-matched");
legend boxoff;
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_Pareto_front.fig",condition_label);
pdf_name = sprintf("%s_Pareto_front.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

%% Plot RMSE and correlation between K and CAQ against flux efficiency
figure;
rmse_K_caq = rmse(log10(K(1:end-1)),log10CAQ_Pareto,1);
plot(log10_flux_scale,rmse_K_caq);
xlabel("log10(\eta/\eta_{max})");
ylabel("RMSE between K and CAQ");
title(condition_label);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_rmse_K_CAQ.fig",condition_label);
pdf_name = sprintf("%s_rmse_K_CAQ.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

figure;
pearson_R_K_caq = corr(log10(K(1:end-1)),log10CAQ_Pareto);
plot(log10_flux_scale,pearson_R_K_caq);
xlabel("log10(\eta/\eta_{max})");
ylabel("Pearson's R between K and CAQ");
title(condition_label);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_pearson_K_CAQ.fig",condition_label);
pdf_name = sprintf("%s_pearson_K_CAQ.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

figure;
spearman_rho_K_caq = corr(log10(K(1:end-1)),log10CAQ_Pareto,'type','Spearman');
plot(log10_flux_scale,spearman_rho_K_caq);
xlabel("log10(\eta/\eta_{max})");
ylabel("Spearman's \rho between K and CAQ");
title(condition_label);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_spearman_K_CAQ.fig",condition_label);
pdf_name = sprintf("%s_spearman_K_CAQ.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

%% Plot FCC values at the CAQ-matched Pareto solution
figure;
bar(FCC_varying_f(:,idx_min_dev_caq));
set(gca,'yscale','log');
axis square;
xticklabels(rxn_names);
title("FCCs, CAQ-matched Pareto solution")
ylabel("FCC");
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_FCC_CAQ_matched.fig",condition_label);
pdf_name = sprintf("%s_FCC_CAQ_matched.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

%% Plot FCC values at the CAQ-matched Pareto solution
figure;
bar(fcc_max_flux);
set(gca,'yscale','log');
axis square;
xticklabels(rxn_names);
title("FCCs, maximal flux")
ylabel("FCC");
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_FCC_max_flux.fig",condition_label);
pdf_name = sprintf("%s_FCC_max_flux.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

%% Plot relationship between optimal E and K
figure;
scatter(K,e_max_flux,64,'Marker','diamond','MarkerEdgeColor',[233 113 50]/256, ...
    'MarkerFaceColor',[233 113 50]/256);
automatic_label(K, e_max_flux, rxn_names, 4, [1.5 0.02]);
xrange = [min(K) max(K)];
yrange = [min(e_max_flux) max(e_max_flux)];
xrange = 10.^[fix(log10(xrange(1)))-1 fix(log10(xrange(2)))+1];
yrange = 10.^[fix(log10(yrange(1)))-1 fix(log10(yrange(2)))+1];
set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel("K");
ylabel("[E], maximal flux");
xlim(xrange);
ylim(yrange);
title(condition_label);
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_K_e_maxflux.fig",condition_label);
pdf_name = sprintf("%s_K_e_maxflux.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);

figure;
scatter(K,e_max_std_varying_f(:,idx_min_dev_caq),64,'Marker','diamond','MarkerEdgeColor',[233 113 50]/256, ...
    'MarkerFaceColor',[233 113 50]/256);
automatic_label(K, e_max_std_varying_f(:,idx_min_dev_caq), rxn_names, 4, [1.5 0.02]);
xrange = [min(K) max(K)];
yrange = [min(e_max_std_varying_f(:,idx_min_dev_caq)) max(e_max_std_varying_f(:,idx_min_dev_caq))];
xrange = 10.^[fix(log10(xrange(1)))-1 fix(log10(xrange(2)))+1];
yrange = 10.^[fix(log10(yrange(1)))-1 fix(log10(yrange(2)))+1];
set(gca,'xscale','log');
set(gca,'yscale','log');
xlabel("K");
ylabel("[E], CAQ-matched");
xlim(xrange);
ylim(yrange);
title(condition_label);
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
fig_name = sprintf("%s_K_e_CAQ_matched.fig",condition_label);
pdf_name = sprintf("%s_K_e_CAQ_matched.pdf",condition_label);
exportgraphics(gcf,pdf_name,'ContentType','vector');
savefig(gcf,fig_name);


%{
%% Plot relationship between rmse of CAQ/[S] and J/J_max
rmse_caq_varying_flux = rmse(log10CAQ_exp(:,1),log10CAQ_Pareto,1);
figure;
plot(log10_flux_scale,rmse_caq_varying_flux);
hold on;
scatter(log10_flux_scale(idx_min_dev_caq),rmse_caq_varying_flux(idx_min_dev_caq), 96,...
    'Marker','diamond','MarkerEdgeColor',[233 113 50]/256,'MarkerFaceColor',[233 113 50]/256);
legend({'Pareto front','CAQ-matched'});
xlabel("log10(J/J_{max})")
ylabel("RMSE for CAQ");
axis square;

rmse_S_varying_flux = rmse(log10(met_conc_exp),log10(S_std_varying_f),1);
figure;
plot(log10_flux_scale,rmse_S_varying_flux);
hold on;
scatter(log10_flux_scale(idx_min_dev_caq),rmse_S_varying_flux(idx_min_dev_caq), 96,...
    'Marker','diamond','MarkerEdgeColor',[233 113 50]/256,'MarkerFaceColor',[233 113 50]/256);
legend({'Pareto front','CAQ-matched'});
xlabel("log10(J/J_{max})")
ylabel("RMSE for [S]");
axis square;

%% Plot the optimization flags
figure;
subplot(1,3,1);
plot(log10_flux_scale,exitflag_std);title("Exitflag");xlabel("log10(J/J_{max})");
subplot(1,3,2);
plot(log10_flux_scale,firstorder_std);
hold on;
scatter(log10_flux_scale(idx_min_dev_caq),firstorder_std(idx_min_dev_caq), 96,...
    'Marker','diamond','MarkerEdgeColor',[233 113 50]/256,'MarkerFaceColor',[233 113 50]/256);
title("Firstorder");xlabel("log10(J/J_{max})");
set(gca,'yscale','log');
subplot(1,3,3);
plot(log10_flux_scale,constrviol_std);
title("Constrviol");xlabel("log10(J/J_{max})");
set(gca,'yscale','log');
%}

%% Definition of nonlinear constraint on flux bound
function [c, ceq] = flux_bound(a,e,K,Sin,Sout,f_min)
    ceq = 0;
    [~,f,~,~]=MCA_Linear(a,e,K,Sin,Sout);
    c = f_min - f; % f>=f_min
end

function [c, ceq] = flux_bound_log10e(a,log10e,K,Sin,Sout,f_min)
    e = 10.^log10e;
    ceq = sum(e) - 1;
    [~,f,~,~]=MCA_Linear(a,e,K,Sin,Sout);
    c = f_min - f; % f>=f_min
end

function s1 = adjust_met_conc(s0,pathway)
% Compute real metabolite concentration from the one simulated from a
% simple linear pathway model that only includes the FBP->DHAP->GAP route
s1 = s0;
if strcmp(pathway,"glycolysis")
    s1(6:end,:) = 2*s0(6:end,:);
elseif strcmp(pathway,"glycolysis_with_GLUT")
    s1(7:end,:) = 2*s0(7:end,:);
end
end

function automatic_label(x_all,y_all,labels,x_dis_cutoff,x_offset)
% automatically optimize label locations to avoid overlapping labels
n_points = length(x_all);
dx = ones(n_points,1)*1.5;
dy = ones(n_points,1);

% adjust locations of labels too close with each other
[~,ord] = sort(log10(y_all),'descend');
for k = 2:length(ord)
    i = ord(k);
    j = ord(k-1);

    if abs(log10(y_all(i)) - log10(y_all(j))) < 0.18 && ...
       abs(log10(x_all(i)) - log10(x_all(j))) < x_dis_cutoff
        if x_all(i) < x_all(j)
            dx(i) = min(x_offset);
            dx(j) = max(x_offset);
        else
            dx(i) = max(x_offset);
            dx(j) = min(x_offset);
        end
    end
end

for i = 1:n_points
    xt = x_all(i) * dx(i);
    yt = y_all(i) * dy(i);

    text(xt, yt, labels{i}, ...
        'FontSize', 20, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
end
end