% ANALYSIS CODE: Gene Expression parameter var and thereshold var

Sin = 1000e3;
%p = 1;
p_range = [1, 3];		% q
%t_off_range = [9, 19, 37, 74];
%t_on_range = [3 6 12];
%threshold_range = [5e5, 1e6, 5e6, 1e7, 2e7];
t_off_range = [4.5, 9, 19, 37];		% j
t_on_range = [6, 12];				% l
threshold_range = [5e6, 1e7, 2e7, 5e7, 1e8];	% i

seed_range = 0:2;

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;
createRec = 0;

same_thresh = 1;

base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_p_exp_thres\sizer\';
%chdir(strcat(base_dir,'\var_exp_threshold\sizer'));

% ADDER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent


%===============
% Extract data from simulation output, and store in usable form
data = {};
data_compiled = {};

for q = 1:length(p_range)
	for i = 1:length(threshold_range)
		for j = 1:length(t_off_range)
			for l = 1:length(t_on_range)					
				data_compiled{q,i,j,l} = struct('div_durs_exp',[],'size_bir_exp',[],'size_div_exp',[]);
				for k = seed_range
					%data{i,j,l,k+1} = load(strcat(base_dir,'var_thr',num2str(threshold_range(i)),'_ton',num2str(t_on_range(l)),'_toff',num2str(t_off_range(j)),'_rng',num2str(k),'.mat'),'div_durs_exp','sim_vars','size_bir','size_div');
					data{q,i,j,l,k+1} = load(strcat(base_dir,'var_thr',num2str(threshold_range(i)),'_ton',num2str(t_on_range(l)),'_toff',num2str(t_off_range(j)),'_p',num2str(p_range(q)),'_rng',num2str(k),'.mat'),'div_durs_exp','sim_vars','config','size_bir_exp','size_div_exp');
					
					data_compiled{q,i,j,l}.div_durs_exp = [data_compiled{q,i,j,l}.div_durs_exp, data{q,i,j,l,k+1}.div_durs_exp(2:end)];
					data_compiled{q,i,j,l}.size_bir_exp = [data_compiled{q,i,j,l}.size_bir_exp; data{q,i,j,l,k+1}.size_bir_exp(2:end,:)];
					data_compiled{q,i,j,l}.size_div_exp = [data_compiled{q,i,j,l}.size_div_exp; data{q,i,j,l,k+1}.size_div_exp(2:end,:)];
				end
			end
		end
	end
end
%save('geneExp_threshold_var_data.mat','data','data_compiled','threshold_range','t_off_range','t_on_range','seed_range','Sin','p');
save('genEx_thres_p_var_data.mat','data_compiled','threshold_range','t_off_range','t_on_range','seed_range','Sin','p_range');%,'data');

%===============
% Load data and plot the effect of the variation

%load('geneExp_threshold_var_data.mat');
load('genEx_thres_p_var_data.mat');

interdiv_stats = zeros(4,length(p_range),length(threshold_range),length(t_off_range),length(t_on_range));
doubling_stats = zeros(4,length(p_range),length(threshold_range),length(t_off_range),length(t_on_range));
doubling_std_stats = zeros(4,length(p_range),length(threshold_range),length(t_off_range),length(t_on_range));

for q = 1:length(p_range)
	for i = 1:length(threshold_range)
		for j = 1:length(t_off_range)
			for l = 1:length(t_on_range)
				interdiv_stats(1,q,i,j,l) = mean(data_compiled{q,i,j,l}.div_durs_exp);
				interdiv_stats(2,q,i,j,l) = std(data_compiled{q,i,j,l}.div_durs_exp);
				interdiv_stats(3,q,i,j,l) = skewness(data_compiled{q,i,j,l}.div_durs_exp);
				interdiv_stats(4,q,i,j,l) = kurtosis(data_compiled{q,i,j,l}.div_durs_exp);
				
				%doubDat = data_compiled{q,i,j,l}.div_durs_exp'./log2(mean(data_compiled{q,i,j,l}.size_div_exp,2)./mean(data_compiled{q,i,j,l}.size_bir_exp,2));
				doubDat = data_compiled{q,i,j,l}.div_durs_exp./log2(mean(data_compiled{q,i,j,l}.size_div_exp')./mean(data_compiled{q,i,j,l}.size_bir_exp'));
				doubling_stats(1,q,i,j,l) = mean(doubDat);
				doubling_stats(2,q,i,j,l) = std(doubDat);
				doubling_stats(3,q,i,j,l) = skewness(doubDat);
				doubling_stats(4,q,i,j,l) = kurtosis(doubDat);
				
				doubDat_std = (doubDat-doubling_stats(1,q,i,j,l))/doubling_stats(2,q,i,j,l);
				doubling_std_stats(1,q,i,j,l) = mean(doubDat_std);
				doubling_std_stats(2,q,i,j,l) = std(doubDat_std);
				doubling_std_stats(3,q,i,j,l) = skewness(doubDat_std);
				doubling_std_stats(4,q,i,j,l) = kurtosis(doubDat_std);
			end
		end
	end
end

%===============
% Variation of threshold

% Plot standardized histograms
q= 2; j = 3; l = 1; i = 3;
figure; hold on;
%for i = length(threshold_range):-1:1
for i = 1:length(threshold_range)
%	for j = 1:length(t_off_range)
%		for l = 1:length(t_on_range)
			%figure; hold on;
			doubDat = data_compiled{q,i,j,l}.div_durs_exp./log2(mean(data_compiled{q,i,j,l}.size_div_exp')./mean(data_compiled{q,i,j,l}.size_bir_exp'));
			histogram((doubDat-doubling_stats(1,q,i,j,l))/doubling_stats(2,q,i,j,l),'normalization','pdf','BinWidth',0.25,'Displayname',strcat('threshold=',num2str(threshold_range(i),'%.0e')));
			pd{i} = fitdist(((doubDat-doubling_stats(1,q,i,j,l))/doubling_stats(2,q,i,j,l))','GeneralizedExtremeValue');
			xl = xlim; yl = ylim;
			%xval = xl(1):diff(xl)/1000:xl(2); 
			xval = xl(1):0.05:12; 
			yval = pdf(pd{i},xval);
			hold on;
			plot(xval,yval,'LineWidth',2,'Displayname','FrechetDist Fit');
			%set(gca, 'YScale', 'log')
%		end
%	end
end
xlabel('Doubling time (in mins)')
ylabel('PDF')
title({'Histogram of standardized Doubling times when Threshold value is altered',''})
legend('location','northeast')
%print('sim_stdDist_varThresh','-r300','-dpng')
% Plot Y in Log scale
set(gca, 'YScale', 'log')
ylim([1e-4 1e0])
print('sim_stdDist_varThresh_LOG','-r300','-dpng')

%===============
% Plot histograms

% Variation of only the threshold
q= 2; j = 3; l = 1; i = 3;
figure; hold on;
%for i = length(threshold_range):-1:1
for i = 1:length(threshold_range)
%	for j = 1:length(t_off_range)
%		for l = 1:length(t_on_range)
			%figure; hold on;
			doubDat = data_compiled{q,i,j,l}.div_durs_exp./log2(mean(data_compiled{q,i,j,l}.size_div_exp')./mean(data_compiled{q,i,j,l}.size_bir_exp'))/60;
			histogram(doubDat,'normalization','pdf','BinWidth',2,'Displayname',strcat('threshold=',num2str(threshold_range(i),'%.0e')));
			pd{i} = fitdist(doubDat','GeneralizedExtremeValue');
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
title({'Histogram of Doubling times varying metab threshold, p = 3',''})
legend('location','northeast')
xlim([0 200])
print('sim_hist_varThresh_p3','-r300','-dpng')
print('sim_hist_varThresh','-r300','-dpng')
% Plot Y in Log scale
set(gca, 'YScale', 'log')
ylim([1e-4 1e-1])
print('sim_hist_varThresh_LOG','-r300','-dpng')

%==========================
% Variation of only t_ON
% Plot standardized histograms
% j = 3; l = 2; i = 5; % tOFF = 37mins, tON var, Sin = 2e7
q= 1; j = 4; l = 1; i = 3;
figure; hold on;
%for i = length(threshold_range):-1:1
%for i = 1:length(threshold_range)
%	for j = 1:length(t_off_range)
		%for l = 1:length(t_on_range)
		for l = length(t_on_range):-1:1
			%figure; hold on;
			doubDat = data_compiled{q,i,j,l}.div_durs_exp./log2(mean(data_compiled{q,i,j,l}.size_div_exp')./mean(data_compiled{q,i,j,l}.size_bir_exp'));
			doubDat_std = (doubDat-doubling_stats(1,q,i,j,l))/doubling_stats(2,q,i,j,l);
			histogram(doubDat_std,'normalization','pdf','BinWidth',0.25,'Displayname',strcat('t_{ON}=',num2str(t_on_range(l),'%.0e')));
			pd{i} = fitdist(doubDat_std','GeneralizedExtremeValue');
			xl = xlim; yl = ylim;
			%xval = xl(1):diff(xl)/1000:xl(2); 
			xval = -2:0.05:12; 
			yval = pdf(pd{i},xval);
			hold on;
			plot(xval,yval,'LineWidth',2,'Displayname','FrechetDist Fit');
			%set(gca, 'YScale', 'log')
%		end
%	end
end
xlabel('Standardized Scale')
ylabel('PDF')
title({'Histogram of standardized Doubling times when t_ON is altered',''})
legend('location','northeast')
print('sim_stdDist_vartON','-r300','-dpng')
% Plot Y in Log scale
set(gca, 'YScale', 'log')
ylim([1e-4 1e-1])
print('sim_stdDist_vartON_LOG','-r300','-dpng')

%==========================
% Variation of only t_ON
% Plot original histograms
%j = 3; l = 2; i = 5; % tOFF = 37mins, tON var, Sin = 2e7
q= 1; j = 4; l = 1; i = 3;
figure; hold on;
%for i = length(threshold_range):-1:1
%for i = 1:length(threshold_range)
%	for j = 1:length(t_off_range)
		%for l = 1:length(t_on_range)
		for l = length(t_on_range):-1:1
			%figure; hold on;
			doubDat = data_compiled{q,i,j,l}.div_durs_exp./log2(mean(data_compiled{q,i,j,l}.size_div_exp')./mean(data_compiled{q,i,j,l}.size_bir_exp'))/60;
			histogram(doubDat,'normalization','pdf','Displayname',strcat('t_{ON}=',num2str(t_on_range(l))));
			pd{i} = fitdist(doubDat','GeneralizedExtremeValue');
			xl = xlim; yl = ylim;
			xval = xl(1):diff(xl)/1000:xl(2); 
			%xval = -2:0.05:12; 
			yval = pdf(pd{i},xval);
			hold on;
			plot(xval,yval,'LineWidth',2,'Displayname','FrechetDist Fit');
			%set(gca, 'YScale', 'log')
%		end
%	end
end
xlabel('Doubling time (in mins)')
ylabel('PDF')
title({'Histogram of Doubling times when t_{ON} is altered',''})
legend('location','northeast')
xlim([0 200])
print('sim_Dist_vartON','-r300','-dpng')
% Plot Y in Log scale
set(gca, 'YScale', 'log')
ylim([1e-4 1e-1])
print('sim_DistLOG_vartON','-r300','-dpng')
%===============

% Variation of only t_OFF
% Plot original histograms
%j = 3; l = 2; i = 5; % tOFF = 37mins, tON var, Sin = 2e7
q= 1; j = 4; l = 1; i = 3;
figure; hold on;
%for i = length(threshold_range):-1:1
%for i = 1:length(threshold_range)
	for j = 1:length(t_off_range)
		%for l = 1:length(t_on_range)
		%for l = length(t_on_range):-1:1
			%figure; hold on;
			doubDat = data_compiled{q,i,j,l}.div_durs_exp./log2(mean(data_compiled{q,i,j,l}.size_div_exp')./mean(data_compiled{q,i,j,l}.size_bir_exp'));
			histogram(doubDat,'normalization','pdf','BinWidth',2,'Displayname',strcat('t_{OFF}= ',num2str(t_off_range(j))));
			pd{i} = fitdist(doubDat','GeneralizedExtremeValue');
			xl = xlim; yl = ylim;
			xval = 0:200/1000:200; 
			%xval = -2:0.05:12; 
			yval = pdf(pd{i},xval);
			hold on;
			plot(xval,yval,'LineWidth',2,'Displayname','FrechetDist Fit');
			%set(gca, 'YScale', 'log')
%		end
%	end
end
xlabel('Doubling time (in mins)')
ylabel('PDF')
title({'Histogram of Doubling times when t_{OFF} is altered',''})
legend('location','northeast')
xlim([0 200])
print('sim_Dist_vartOFF','-r300','-dpng')
% Plot Y in Log scale
set(gca, 'YScale', 'log')
ylim([1e-4 1e-1])
print('sim_DistLOG_vartOFF','-r300','-dpng')

%===============

% Variation of only t_OFF
% Plot standardized histograms
%j = 3; l = 2; i = 5; % tOFF = 37mins, tON var, Sin = 2e7
q= 1; j = 4; l = 1; i = 3;
figure; hold on;
%for i = length(threshold_range):-1:1
%for i = 1:length(threshold_range)
	%for j = 1:length(t_off_range)
	for j = length(t_off_range):-1:1
		%for l = 1:length(t_on_range)
		%for l = length(t_on_range):-1:1
			%figure; hold on;
			doubDat = data_compiled{q,i,j,l}.div_durs_exp./log2(mean(data_compiled{q,i,j,l}.size_div_exp')./mean(data_compiled{q,i,j,l}.size_bir_exp'));
			doubDat_std = (doubDat-doubling_stats(1,q,i,j,l))/doubling_stats(2,q,i,j,l);
			histogram(doubDat_std,'normalization','pdf','Displayname',strcat('t_{OFF}= ',num2str(t_off_range(j))));
			pd{i} = fitdist(doubDat_std','GeneralizedExtremeValue');
			xl = xlim; yl = ylim;
			xval = xl(1):diff(xl)/1000:xl(2); 
			%xval = -2:0.05:12; 
			yval = pdf(pd{i},xval);
			hold on;
			plot(xval,yval,'LineWidth',2,'Displayname','FrechetDist Fit');
			%set(gca, 'YScale', 'log')
%		end
%	end
end
xlabel('Standardized Scale')
ylabel('PDF')
title({'Histogram of Doubling times when t_{OFF} is altered',''})
legend('location','northeast')
xlim([0 200])
print('sim_stdDist_vartOFF','-r300','-dpng')
% Plot Y in Log scale
set(gca, 'YScale', 'log')
ylim([1e-4 1e0])
print('sim_stdDistLOG_vartOFF','-r300','-dpng')










% ======================
% Plot change in statistics of histogram

% Variation of threshold, at different gene expression parameters
j = 3; % t_off 3=> 37 min
l = 2; % t_ON 2=> 6 min
figure; hold on;
for j = 1:length(t_off_range)
%for l = 1:length(t_on_range)
	%plot(threshold_range, doubling_stats(1,:,j,l)'/60,'-o','linewidth',2);
	errorbar(threshold_range, doubling_stats(1,:,j,l)'/60, doubling_stats(2,:,j,l)'/60,'-o','linewidth',2,'MarkerSize',10);
end
%end
xlabel('Increasing metabolite threshold')
ylabel('Mean doubling time (in mins)')

title('Variation in the Moments of the distributions')
legend('Mean','Std Dev','Skewness','Kurtosis','location','northwest')


% Variation of only the threshold, and standardized data
j = 3; l = 2;
figure; plot(doubling_std_stats([1 2 3 4],:,j,l)','-o','linewidth',2);
xlabel('Increasing metabolite threshold')
ylabel('Value')
title('Variation in the Moments of the distributions')
legend('Mean','Std Dev','Skewness','Kurtosis','location','best')

