nrxn=10;
nmet=nrxn-1;
nsamples=10000; %Number of sampled parameter sets

% Read data from processed realistic datasets
Keq_data = readmatrix("data/K.csv");
kcat_data = readmatrix("data/kcat.csv");
Km_data = readmatrix("data/Km.csv");
S_data = readmatrix("data/S.csv");

% Randomly retrieve data points from the realistic datasets
Keq = reshape(Keq_data(randi(numel(Keq_data), nrxn*nsamples, 1)), nrxn, nsamples);
kcat = reshape(kcat_data(randi(numel(kcat_data), nrxn*nsamples, 1)), nrxn, nsamples);
Ks = reshape(Km_data(randi(numel(Km_data), nrxn*nsamples, 1)), nrxn, nsamples);
Kp = reshape(Km_data(randi(numel(Km_data), nrxn*nsamples, 1)), nrxn, nsamples);
Sin = S_data(randi(numel(S_data), 1, nsamples));
Sout = S_data(randi(numel(S_data), 1, nsamples));

%count=0;
c = zeros(nsamples,1);
cadj = c;
rmse_k_caq = c;
rmse_k_caq_adj = c;
J = zeros(nsamples,1);
Jflag = J;
eMat = zeros(nsamples,nrxn);
phiMat = zeros(nsamples,nrxn-1);
ae2rMat = phiMat;
saturationMat = eMat;

%% simulate MM kinetics with randomly sampled parameters
for i = 1:nsamples
    i
    dG_overall=log(Sout(i)/Sin(i)/prod(Keq(:,i)));
    if dG_overall>0
        b=Sin(i);
        Sin(i)=Sout(i);
        Sout(i)=b;
        kcat(:,i)=kcat(nrxn:-1:1,i);
        Ks(:,i)=Ks(nrxn:-1:1,i);
        Kp(:,i)=Kp(nrxn:-1:1,i);
        Keq(:,i)=1./Keq(nrxn:-1:1,i);
    end
    [J(i),e,Jflag(i)] = maxFluxEfficiency(Sin(i),Sout(i),kcat(:,i),Ks(:,i),Kp(:,i),Keq(:,i));    
    [c(i),cadj(i),rmse_k_caq(i),rmse_k_caq_adj(i),sat,phi,ae2r] = corr_Keq_aEratio(Keq(:,i), ...
        kcat(:,i),Ks(:,i),Kp(:,i),Sin(i),Sout(i),e);
    eMat(i,:) = reshape(e,1,nrxn);
    phiMat(i,:) = reshape(phi,1,nrxn-1);
    ae2rMat(i,:) = reshape(ae2r,1,nrxn-1);
    saturationMat(i,:) = reshape(sat,1,nrxn);
end 

%% Compute c and cadj again using the log10 transformed values
for i = 1:nsamples
    c(i) = corr(log10(Keq(1:9,i)),log10(ae2rMat(i,:)'));
    cadj(i) = corr(log10(Keq(1:9,i).*phiMat(i,:)'),log10(ae2rMat(i,:)'));
end

figure;
histogram(c(J>0 & Jflag>0), 20, 'LineWidth',1,'EdgeColor','k');
xlabel("Pearson's R");
ylabel("Number of models");
title(sprintf("Distribution of Pearson's R\n between log10(K) and log10(CAQ)"));
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
exportgraphics(gcf,'dist_Pearson_R.pdf','ContentType','vector');
savefig(gcf,'dist_Pearson_R.fig');

%% Plot histogram of RMSE
figure;
histogram(rmse_k_caq(J>0 & Jflag>0), 20, 'LineWidth',1,'EdgeColor','k');
xlabel("RMSE");
ylabel("Number of models");
title(sprintf("Distribution of RMSE\n between log10(K) and log10(CAQ)"));
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
exportgraphics(gcf,'dist_logRMSE.pdf','ContentType','vector');
savefig(gcf,'dist_logRMSE.fig');

%% Analyze the relationship between the K-aE2 correlation and saturation terms
cmap = brewermap(256,'RdYlBu');
cmap = cmap(end:-1:1,:);
figure;
dscatter(mean(saturationMat(J>0 & Jflag>0,:),2),c(J>0 & Jflag>0));
xlabel("Mean saturation term");
ylabel("Pearson's R");
box on;
title("Pearson's R between log10(K) and log10(CAQ)");
hold on;
plot([0 1],[0 0],":k");
colormap(cmap);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
exportgraphics(gcf,'saturation_vs_Pearson_R.pdf','ContentType','vector');
savefig(gcf,'saturation_vs_Pearson_R.fig');

%% Analyze the relationship between the K-aE2 RMSE and saturation terms
figure;
dscatter(mean(saturationMat(J>0 & Jflag>0,:),2),rmse_k_caq(J>0 & Jflag>0));
hold on;
plot([0 1],[1 1],":k");
xlabel("Mean saturation term");
ylabel("RMSE");
box on;
title("RMSE between log10(K) and log10(CAQ)");
colormap(cmap);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
%set(gca,'YScale','log');
exportgraphics(gcf,'saturation_vs_RMSE.pdf','ContentType','vector');
savefig(gcf,'saturation_vs_RMSE.fig');

%% Compare the Pearson's R, unadjusted vs adjusted
figure;
%subplot(1,2,1);
dscatter(c(J>0 & Jflag>0),cadj(J>0 & Jflag>0));
hold on;
plot([-1 1],[-1 1],':k');
xlabel("log10(K) vs log10(CAQ)");
ylabel("log10(K\Phi) vs log10(CAQ)");
box on;
title(sprintf("Pearson's R between log10(K) or\n log10(K\\Phi) and log10(CAQ)"));
xlim([0 1]);
ylim([0 1]);
colormap(cmap);
axis square;
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
exportgraphics(gcf,'Pearson_R_KPhi.pdf','ContentType','vector');
savefig(gcf,'Pearson_R_KPhi.fig');

%% Compare the RMSE, unadjusted vs adjusted
figure;
x = (rmse_k_caq(J>0 & Jflag>0));
y = (rmse_k_caq_adj(J>0 & Jflag>0));
range = [min([x;y]) max([x;y])];
dscatter(x,y);
hold on;plot(range,range,":k");
xlabel("log10(K) vs log10(CAQ)");
ylabel("log10(K\Phi) vs log10(CAQ)");
box on;
title(sprintf("RMSE between log10(K) or\n log10(K\\Phi) and log10(CAQ)"));
axis square;
xlim(range);
ylim(range);
colormap(cmap);
set(gcf,'Units','centimeters','Position',[0 0 21 21]);
set(gca, 'Units','centimeters','Position',[3 3 13 13]);
exportgraphics(gcf,'logRMSE_KPhi.pdf','ContentType','vector');
savefig(gcf,'logRMSE_KPhi.fig');

%% Plot some K vs CAQ
figure;
idx = randi(length(idx_good),25);
for i = 1:25
    subplot(5,5,i);
    scatter(log10(Keq(1:9,idx(i))),log10(ae2rMat(idx(i),:)'));
    title(sprintf("c=%.2f,rmse=%.2f",c(idx(i)),rmse_k_caq(idx(i))));
    xlabel("log10K");
    ylabel("log10CAQ");
end

