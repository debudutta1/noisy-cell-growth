% Gene Expression parameter var and thereshold var

Sin = 10000e3;
p = 1;

t_off_range = 37;
t_on_range = 6;

%t_off_range = [19 74];
%t_on_range = [3 12];

threshold_range = [1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9];

seed_range = 0:2;

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;
createRec = 0;

same_thresh = 1;

base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';

%===============
% Extract data from simulation output, and store in usable form
% Sizer
data_sizer = {};
data_compiled_sizer = {};

for i = 1:length(threshold_range)
	for j = 1:length(t_off_range)
		for l = 1:length(t_on_range)					
			data_compiled_sizer{i} = struct('div_durs_exp',[],'size_bir',[],'size_div',[]);
			for k = seed_range
				data_sizer{i,k+1} = load(strcat(base_dir,'\var_threshold\sizer\','var_thr',num2str(threshold_range(i)),'_ton',num2str(t_on_range(l)),'_toff',num2str(t_off_range(j)),'_rng',num2str(k),'.mat'),'div_durs_exp','sim_vars','size_bir','size_div');
				data_compiled_sizer{i}.div_durs_exp = [data_compiled_sizer{i}.div_durs_exp, data_sizer{i,k+1}.div_durs_exp];
				data_compiled_sizer{i}.size_bir = [data_compiled_sizer{i}.size_bir; data_sizer{i,k+1}.size_bir];
				data_compiled_sizer{i}.size_div = [data_compiled_sizer{i}.size_div; data_sizer{i,k+1}.size_div];
			end
		end
	end
end

%Adder
data_adder = {};
data_compiled_adder = {};

for i = 1:length(threshold_range)
	for j = 1:length(t_off_range)
		for l = 1:length(t_on_range)					
			data_compiled_adder{i} = struct('div_durs_exp',[],'size_bir',[],'size_div',[]);
			for k = seed_range
				data_adder{i,k+1} = load(strcat(base_dir,'\var_threshold\adder\','var_thr',num2str(threshold_range(i)),'_ton',num2str(t_on_range(l)),'_toff',num2str(t_off_range(j)),'_rng',num2str(k),'.mat'),'div_durs_exp','sim_vars','size_bir','size_div');
				data_compiled_adder{i}.div_durs_exp = [data_compiled_adder{i}.div_durs_exp, data_adder{i,k+1}.div_durs_exp];
				data_compiled_adder{i}.size_bir = [data_compiled_adder{i}.size_bir; data_adder{i,k+1}.size_bir];
				data_compiled_adder{i}.size_div = [data_compiled_adder{i}.size_div; data_adder{i,k+1}.size_div];
			end
		end
	end
end
save('threshold_var_data.mat','data_sizer','data_compiled_sizer','data_adder','data_compiled_adder','threshold_range','t_off_range','t_on_range','seed_range','Sin','p');

%===============================

% Load data and plot the effect of the variation
load('threshold_var_data.mat');

% Sizer
sz_interdiv_stats = zeros(4,length(threshold_range));
sz_doubling_stats = zeros(4,length(threshold_range));
sz_doubling_std_stats = zeros(4,length(threshold_range));
sz_doubling_norm_stats = zeros(4,length(threshold_range));

for i = 1:length(threshold_range)
	sz_interdiv_stats(1,i) = mean(data_compiled_sizer{i}.div_durs_exp);
	sz_interdiv_stats(2,i) = std(data_compiled_sizer{i}.div_durs_exp);
	sz_interdiv_stats(3,i) = skewness(data_compiled_sizer{i}.div_durs_exp);
	sz_interdiv_stats(4,i) = kurtosis(data_compiled_sizer{i}.div_durs_exp);
			
	doubDat = data_compiled_sizer{i}.div_durs_exp./log2(mean(data_compiled_sizer{i}.size_div)./mean(data_compiled_sizer{i}.size_bir));
	sz_doubling_stats(1,i) = mean(doubDat);
	sz_doubling_stats(2,i) = std(doubDat);
	sz_doubling_stats(3,i) = skewness(doubDat);
	sz_doubling_stats(4,i) = kurtosis(doubDat);
			
	doubDat_std = (doubDat-sz_doubling_stats(1,i))/sz_doubling_stats(2,i);
	sz_doubling_std_stats(1,i) = mean(doubDat_std);
	sz_doubling_std_stats(2,i) = std(doubDat_std);
	sz_doubling_std_stats(3,i) = skewness(doubDat_std);
	sz_doubling_std_stats(4,i) = kurtosis(doubDat_std);
	
	doubDat_norm = doubDat/sz_doubling_stats(1,i);
	sz_doubling_norm_stats(1,i) = mean(doubDat_norm);
	sz_doubling_norm_stats(2,i) = std(doubDat_norm);
	sz_doubling_norm_stats(3,i) = skewness(doubDat_norm);
	sz_doubling_norm_stats(4,i) = kurtosis(doubDat_norm);
end
% Adder
ad_interdiv_stats = zeros(4,length(threshold_range));
ad_doubling_stats = zeros(4,length(threshold_range));
ad_doubling_std_stats = zeros(4,length(threshold_range));
ad_doubling_norm_stats = zeros(4,length(threshold_range));

for i = 1:length(threshold_range)
	ad_interdiv_stats(1,i) = mean(data_compiled_adder{i}.div_durs_exp);
	ad_interdiv_stats(2,i) = std(data_compiled_adder{i}.div_durs_exp);
	ad_interdiv_stats(3,i) = skewness(data_compiled_adder{i}.div_durs_exp);
	ad_interdiv_stats(4,i) = kurtosis(data_compiled_adder{i}.div_durs_exp);
			
	doubDat = data_compiled_adder{i}.div_durs_exp./log2(mean(data_compiled_adder{i}.size_div)./mean(data_compiled_adder{i}.size_bir));
	ad_doubling_stats(1,i) = mean(doubDat);
	ad_doubling_stats(2,i) = std(doubDat);
	ad_doubling_stats(3,i) = skewness(doubDat);
	ad_doubling_stats(4,i) = kurtosis(doubDat);
			
	doubDat_std = (doubDat-ad_doubling_stats(1,i))/ad_doubling_stats(2,i);
	ad_doubling_std_stats(1,i) = mean(doubDat_std);
	ad_doubling_std_stats(2,i) = std(doubDat_std);
	ad_doubling_std_stats(3,i) = skewness(doubDat_std);
	ad_doubling_std_stats(4,i) = kurtosis(doubDat_std);
	
	doubDat_norm = doubDat/ad_doubling_stats(1,i);
	ad_doubling_norm_stats(1,i) = mean(doubDat_norm);
	ad_doubling_norm_stats(2,i) = std(doubDat_norm);
	ad_doubling_norm_stats(3,i) = skewness(doubDat_norm);
	ad_doubling_norm_stats(4,i) = kurtosis(doubDat_norm);
end

%===============================

figure; hold on;
%plot(threshold_range, sz_doubling_stats(1,:)'/60,'-o','linewidth',2);
%plot(threshold_range, ad_doubling_stats(1,:)'/60,'-o','linewidth',2);
errorbar(threshold_range, sz_doubling_stats(1,:)'/60, sz_doubling_stats(2,:)'/60,'-o','linewidth',2)%,'MarkerSize',10);
errorbar(threshold_range, ad_doubling_stats(1,:)'/60, ad_doubling_stats(2,:)'/60,'-o','linewidth',2)%,'MarkerSize',10);
xlabel('Increasing metabolite threshold')
ylabel('Mean doubling time (in mins)')


figure; hold on;
%plot(threshold_range, sz_doubling_stats(1,:)'/60,'-o','linewidth',2);
%plot(threshold_range, ad_doubling_stats(1,:)'/60,'-o','linewidth',2);
plot(threshold_range, sz_doubling_stats([1:4],:)'/60,'-o','linewidth',2)%,'MarkerSize',10);
%plot(threshold_range, ad_doubling_norm_stats(2,:)'/60,'-o','linewidth',2)%,'MarkerSize',10);
xlabel('Increasing metabolite threshold')
ylabel('Mean doubling time (in mins)')

%===============================
% Plot histograms standardized

figure; hold on;
%for i = 1:5
%for i = length(threshold_range):-1:1
for i = 1:length(threshold_range)
			%figure;
			doubDat = data_compiled_sizer{i}.div_durs_exp./log2(mean(data_compiled_sizer{i}.size_div)./mean(data_compiled_sizer{i}.size_bir));
			doubDat_std = (doubDat-ad_doubling_stats(1,i))/ad_doubling_stats(2,i);
			
			% Original data
			histogram(doubDat/60,'normalization','pdf','Displayname',strcat('threshold=',num2str(threshold_range(i),'%.0e')));
			pd{i} = fitdist(doubDat','GeneralizedExtremeValue');
			
			% standardized data
%			histogram(doubDat_std,'normalization','pdf','BinWidth',0.25,'Displayname',strcat('threshold=',num2str(threshold_range(i),'%.0e')));
			%histogram(doubDat_std,'normalization','pdf','Displayname',strcat('threshold=',num2str(threshold_range(i),'%.0e')));
%			pd{i} = fitdist(doubDat_std','GeneralizedExtremeValue');
			
			% Plot FrechetDist
			xl = xlim; yl = ylim;
			xval = xl(1):diff(xl)/1000:xl(2); 
			%xval = -2:0.05:12; 
			yval = pdf(pd{i},xval);
			hold on;
			plot(xval,yval,'LineWidth',2,'Displayname','FrechetDist Fit');
			%set(gca, 'YScale', 'log')
			%xlabel('Doubling time (in mins)')
			%ylabel('PDF')
%		end
%	end
end
xlabel('Doubling time (in mins)')
ylabel('PDF')
title({'Histogram of Doubling times when Threshold value is altered',''})
legend('location','northeast')

print('sim_stdDst_std_varThresh','-r300','-dpng')
%print('sim_hist_varThresh','-r300','-dpng')

%=================
% Plot 3D Histograms stacked using Waterfall plot
Z = [];
Y = 1:length(threshold_range); % Number of data points and how separated we want them to be
hist_edges = -3.2:0.1:8.4;
%X = % Mid point of bin edges
tmp = movmean(hist_edges,2);
X = tmp(2:end);

for i = 1:length(threshold_range)
	doubDat = data_compiled_sizer{i}.div_durs_exp./log2(mean(data_compiled_sizer{i}.size_div)./mean(data_compiled_sizer{i}.size_bir));
	doubDat_std = (doubDat-ad_doubling_stats(1,i))/ad_doubling_stats(2,i);
	[Z(i,:)] = histcounts(doubDat_std,hist_edges,'normalization','pdf');
end
figure;
% Waterfall
ax1 = waterfall(X,Y,Z);
set(ax1, 'Linewidth', 3)
% Ribbon
figure;
ax2 = ribbon(X,Z');

%=================
%=================
% Plot the mean doubling time with threshold
figure; hold on;
plot(threshold_range, sz_doubling_stats(1,:)'/60,'-o','linewidth',2);
%errorbar(threshold_range, sz_doubling_stats(1,:)'/60, sz_doubling_stats(2,:)'/60,'-o','linewidth',2,'MarkerSize',10);
set(gca, 'XScale', 'log')
xlabel('Increasing metabolite threshold')
ylabel('Mean doubling time (in mins)')

% Plot the Frechet Fit parameters for each of the standardized distributions for varying metabolite threshold
freschet_sz = struct('k',[],'sig',[],'mu',[]);
for i = 1:length(threshold_range)
	freschet_sz.k(i) = pd{i}.k;
	freschet_sz.sig(i) = pd{i}.sigma;
	freschet_sz.mu(i) = pd{i}.mu;
end
figure; plot(threshold_range, freschet_sz.k,'-o','linewidth',2); 
set(gca, 'XScale', 'log');
xlabel('Increasing metabolite threshold')
ylabel('Frechet fit shape parameter (k)')

figure; plot(threshold_range, freschet_sz.sig,'-o','linewidth',2); 
set(gca, 'XScale', 'log');
xlabel('Increasing metabolite threshold')
ylabel('Frechet fit scale parameter (sigma)')

figure; plot(threshold_range, freschet_sz.mu,'-o','linewidth',2); 
set(gca, 'XScale', 'log');
xlabel('Increasing metabolite threshold')
ylabel('Frechet fit location parameter (mu)')