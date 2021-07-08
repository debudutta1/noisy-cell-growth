% load the data
p = 3;
n_ext_range = [1:p];
seed_range = 0:4;
%feed_rate = [0:2500:15000];
feed_rate = [0 625 1250 2500:2500:15000];
Sin = 1e3*[7.5 10 13 20 50];

% ==============
%% (i) Shifted Histogram, p=3, high Sin=50k,and low sin= 7.5k
% ==============
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_extFeed\sizer_toff_37\';
load(strcat(base_dir,'var_extFeed_fig4.mat'))

% RESCALED TO MEAN ==== Varying Import

%% Low Sin
n_ext = 1;
f = figure; 
plot_indices = 1:2:9;
subplot(2,1,1); hold on;
i = 1; 
%for j = 1:length(feed_rate)
for j = 1:length(plot_indices) %1:length(feed_rate)
	doubDat = div_durs_compiled(:,:,j,n_ext,i)/60;
	doubDat = doubDat(:);
	doubDat = doubDat/mean(doubDat);
	
	dat_std(j) = std(doubDat);
	dat_skew(j) = skewness(doubDat);
	dat_kurt(j) = kurtosis(doubDat);
	
	pd = fitdist(doubDat(:),'GeneralizedExtremeValue');
	gev_k(j) = pd.k;
	gev_sig(j) = pd.sigma;
	gev_mu(j) = pd.mu;
	
end

f = figure; hold on;
subplot(3,1,1)
plot(feed_rate(plot_indices),dat_std,'-o','linewidth',2)
ylim([0.42 0.5])
ylabel('CV','FontSize',11)
subplot(3,1,2)
plot(feed_rate(plot_indices),dat_skew,'-o','linewidth',2)
ylabel('Skewness','FontSize',11)
subplot(3,1,3)
plot(feed_rate(plot_indices),dat_kurt,'-o','linewidth',2)
ylabel('Kurtosis','FontSize',11)
xlabel('Varying Import Rate (per sec)','FontSize',11)

%sgtitle({'Rescaled histograms of generation time, for varying p, when Sin = 50k/sec',''})

name = 'scaled_Hist_varFeed_lowSin_SUPP_CV_skew';
set(f, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r600', '-transparent', '-painters')

%======================
%% High Sin
n_ext = 1;
f = figure; 
plot_indices = 1:2:9;
subplot(2,1,1); hold on;
i = 5; 
%for j = 1:length(feed_rate)
for j = 1:length(plot_indices) %1:length(feed_rate)
	doubDat = div_durs_compiled(:,:,j,n_ext,i)/60;
	doubDat = doubDat(:);
	doubDat = doubDat/mean(doubDat);
	
	dat_std(j) = std(doubDat);
	dat_skew(j) = skewness(doubDat);
	dat_kurt(j) = kurtosis(doubDat);
	
	pd = fitdist(doubDat(:),'GeneralizedExtremeValue');
	gev_k(j) = pd.k;
	gev_sig(j) = pd.sigma;
	gev_mu(j) = pd.mu;
	
end

f = figure; hold on;
subplot(3,1,1)
plot(feed_rate(plot_indices),dat_std,'-o','linewidth',2)
ylim([0.42 0.46])
ylabel('CV','FontSize',11)
subplot(3,1,2)
plot(feed_rate(plot_indices),dat_skew,'-o','linewidth',2)
ylabel('Skewness','FontSize',11)
subplot(3,1,3)
plot(feed_rate(plot_indices),dat_kurt,'-o','linewidth',2)
ylabel('Kurtosis','FontSize',11)
xlabel('Varying Import Rate (per sec)','FontSize',11)

%sgtitle({'Rescaled histograms of generation time, for varying p, when Sin = 50k/sec',''})

name = 'scaled_Hist_varFeed_highSin_SUPP_CV_skew';
set(f, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r600', '-transparent', '-painters')

%=============================

%% load data
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';
load(strcat(base_dir,'\protSec\Synced\Sizer\','var_protSec_fig4_p3.mat'),'div_durs_compiled','growth_rate');

% RESCALED TO MEAN ==== Varying Secretion
Sin = 1e3*[7.5 10 13 20 50];
p = 3;
seed_range = 0:2;
secRatio = [0 0.05 0.1 0.25 0.5 0.75];


%% Low Sin
i = 1; j = 1; 	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	doubDat = div_durs_compiled(:,:,j,l,i)/60;
	doubDat = doubDat(:);
	doubDat = doubDat/mean(doubDat);
	
	dat_std(l) = std(doubDat);
	dat_skew(l) = skewness(doubDat);
	dat_kurt(l) = kurtosis(doubDat);
	
	pd = fitdist(doubDat(:),'GeneralizedExtremeValue');
	gev_k(l) = pd.k;
	gev_sig(l) = pd.sigma;
	gev_mu(l) = pd.mu;
	
end

f = figure; hold on;
subplot(3,1,1)
plot(secRatio,dat_std,'-o','linewidth',2)
%ylim([0.42 0.5])
ylabel('CV','FontSize',11)
subplot(3,1,2)
plot(secRatio,dat_skew,'-o','linewidth',2)
ylabel('Skewness','FontSize',11)
subplot(3,1,3)
plot(secRatio,dat_kurt,'-o','linewidth',2)
ylabel('Kurtosis','FontSize',11)
xlabel('Varying Secrete Rate (per sec)','FontSize',11)

%sgtitle({'Rescaled histograms of generation time, for varying p, when Sin = 50k/sec',''})

name = 'scaled_Hist_varSecrete_lowSin_SUPP_CV_skew';
set(f, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r600', '-transparent', '-painters')

%===============================

%% High Sin
i = 5; j = 1; 	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	doubDat = div_durs_compiled(:,:,j,l,i)/60;
	doubDat = doubDat(:);
	doubDat = doubDat/mean(doubDat);
	
	dat_std(l) = std(doubDat);
	dat_skew(l) = skewness(doubDat);
	dat_kurt(l) = kurtosis(doubDat);
	
	pd = fitdist(doubDat(:),'GeneralizedExtremeValue');
	gev_k(l) = pd.k;
	gev_sig(l) = pd.sigma;
	gev_mu(l) = pd.mu;
	
end

f = figure; hold on;
subplot(3,1,1)
plot(secRatio,dat_std,'-o','linewidth',2)
%ylim([0.42 0.5])
ylabel('CV','FontSize',11)
subplot(3,1,2)
plot(secRatio,dat_skew,'-o','linewidth',2)
ylabel('Skewness','FontSize',11)
subplot(3,1,3)
plot(secRatio,dat_kurt,'-o','linewidth',2)
ylabel('Kurtosis','FontSize',11)
xlabel('Varying Secrete Rate (per sec)','FontSize',11)

%sgtitle({'Rescaled histograms of generation time, for varying p, when Sin = 50k/sec',''})

name = 'scaled_Hist_varSecrete_highSin_SUPP_CV_skew';
set(f, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r600', '-transparent', '-painters')