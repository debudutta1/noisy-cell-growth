% NEW Analysis Code - 20 March 2020

% New Adder vs Sizer

p = [1 2 3 5 10];
Sin = 1e3*[5 7.5 11 20 50];
seed_range = 0:4;

% Load DATA
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_p_sin\Synced\';

% ADDER
for j = 1:length(p)
	for i = 1:length(Sin)
		for k = seed_range
			clear div_durs div_durs_exp sim_vars size_bir_exp size_div_exp
			load(strcat(base_dir,'Adder\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','sim_vars','size_bir_exp', 'size_div_exp');
			
			x_dat{k+1,i,j} = div_durs_exp;
			
			x_bsiz{k+1,i,j} = size_bir_exp;
			x_dsiz{k+1,i,j} = size_div_exp;
			
			x_mean(k+1,i,j) = mean(div_durs_exp);
			
			x_simvars{k+1,i,j} = sim_vars;
			%x_median(k+1,i,j) = median(div_durs_exp);
			%x_q25(k+1,i,j) = quantile(div_durs_exp,0.25);
			%x_q75(k+1,i,j) = quantile(div_durs_exp,0.75);
		end
	end
end
save('sin_p_var_adderDat.mat')

% SIZER
for j = 1:length(p)
	for i = 1:length(Sin)
		for k = seed_range
			clear div_durs div_durs_exp sim_vars size_bir_exp size_div_exp
			load(strcat(base_dir,'Sizer\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','sim_vars','size_bir_exp', 'size_div_exp');
			
			x_dat{k+1,i,j} = div_durs_exp;
			
			x_bsiz{k+1,i,j} = size_bir_exp;
			x_dsiz{k+1,i,j} = size_div_exp;
			
			x_mean(k+1,i,j) = mean(div_durs_exp);
			
			x_simvars{k+1,i,j} = sim_vars;
			%x_median(k+1,i,j) = median(div_durs_exp);
			%x_q25(k+1,i,j) = quantile(div_durs_exp,0.25);
			%x_q75(k+1,i,j) = quantile(div_durs_exp,0.75);
		end
	end
end
save('sin_p_var_sizerDat.mat')

%========================================
% Load data and process
sizer = load('sin_p_var_adderDat.mat');
adder = load('sin_p_var_sizerDat.mat');

% Create standardized datasets
for j = 1:length(p)
	for i = 1:length(Sin)
		for k = seed_range+1
			sizer.std_dat{k,i,j} = (sizer.x_dat{k,i,j} - mean(sizer.x_dat{k,i,j}))./std(sizer.x_dat{k,i,j});
			tem_dat = sizer.x_dat{k,i,j}'/60./log2(mean(sizer.x_dsiz{k,i,j},2)./mean(sizer.x_bsiz{k,i,j},2));
			tem_dat = tem_dat(2:end);	% Neglect the 1st data point has size=0!
			sizer.std_doub{k,i,j} = (tem_dat - mean(tem_dat))./std(tem_dat);
			
			adder.std_dat{k,i,j} = (adder.x_dat{k,i,j} - mean(adder.x_dat{k,i,j}))./std(adder.x_dat{k,i,j});
			tem_dat = adder.x_dat{k,i,j}'/60./log2(mean(adder.x_dsiz{k,i,j},2)./mean(adder.x_bsiz{k,i,j},2));
			tem_dat = tem_dat(2:end);	% Neglect the 1st data point has size=0!
			adder.std_doub{k,i,j} = (tem_dat - mean(tem_dat))./std(tem_dat);
		end
	end
end

%========================================
% Plot inter-division time histogram
figure;histogram(x_dat{3,2,2}/60,'normalization','pdf')

% Plot Doubling time histogram
figure;histogram(x_dat{3,2,2}/60./log2(mean(x_dsiz{3,2,2},2)./mean(x_bsiz{3,2,2},2))','normalization','pdf')

%========================================
% Compare Adder vs Sizer

% Histogram of inter-division times
figure; hold on;
histogram(sizer.x_dat{3,2,2}/60,'normalization','pdf','DisplayName','Sizer')
histogram(adder.x_dat{3,2,2}/60,'normalization','pdf','DisplayName','Adder')
xlabel('Cell inter-division times (in mins)')
ylabel('PDF')
title({'Histogram of inter-division times, from Adder vs Sizer model based simulations',''})
legend('location','northeast')
print('InterDiv_dist_ADvsSZ','-dpng','-r300')


figure; hold on;
k = 3; i = 2; j = 2;
tem_dat = sizer.x_dat{k,i,j}'/60./log2(mean(sizer.x_dsiz{k,i,j},2)./mean(sizer.x_bsiz{k,i,j},2));
histogram(tem_dat,'normalization','pdf','DisplayName','Sizer')
tem_dat = adder.x_dat{k,i,j}'/60./log2(mean(adder.x_dsiz{k,i,j},2)./mean(adder.x_bsiz{k,i,j},2));
histogram(tem_dat,'normalization','pdf','DisplayName','Adder')
xlabel('Cell doubling times (in mins)')
ylabel('PDF')
title({'Histogram of cell Doubling times, from Adder vs Sizer model based simulations',''})
legend('location','northeast')
print('Doubling_dist_ADvsSZ','-dpng','-r300')


%========================================
%========================================
%========================================

% Plot histogram with varying sin, for each p
k = 1;
for j = 1:length(p)	
	figure; hold on;
	for i = 1:length(Sin)
		%histogram(x_dat{1,i,j}/60,'normalization','pdf');
		%histogram(x_dat{k,i,j}/60./log2(mean(x_dsiz{k,i,j},2)./mean(x_bsiz{k,i,j},2))','normalization','pdf')
		
		%histogram(sizer.x_dat{k,i,j}/60./log2(mean(sizer.x_dsiz{k,i,j},2)./mean(sizer.x_bsiz{k,i,j},2))','normalization','pdf')
		%histogram(sizer.std_doub{k,i,j},'normalization','pdf')
		
		%histogram(adder.x_dat{k,i,j}/60./log2(mean(adder.x_dsiz{k,i,j},2)./mean(adder.x_bsiz{k,i,j},2))','normalization','pdf')
		histogram(adder.std_doub{k,i,j},'normalization','pdf')
	end
end

% Plot histogram with varying p, for each Sin
for i = 1:length(Sin)
	figure; hold on;
	for j = 1:length(p)	
		%histogram(x_dat{1,i,j}/60,'normalization','pdf');
		histogram(x_dat{k,i,j}/60./log2(mean(x_dsiz{k,i,j},2)./mean(x_bsiz{k,i,j},2))','normalization','pdf')
	end
end