% Load data - SUPP for Fig 3b

%% load data
base_dir = 'D:\Debu Simulations\Sep 2020\var_protSec\';
load(strcat(base_dir,'var_protSec_fig4.mat'),'div_durs_compiled','growth_rate');

Sin = 1e3*[1.5 5 6 6.5 7 9 13];
% p = [1 3 5];
p = 3;
seed_range = 0:2;
secRatio = [0 0.05 0.1 0.25 0.5 0.75];

lin_colors = lines(14);

% ========================================
% STANDARDIZED => RESCALED TO MEAN

%i = 2	%low sin
%i = 7	% high sin

j = 1; 	% j = 1 => p = 3. Only p that data extracted

% div_durs_compiled(:,:,j,l,i)		=> l - secRatio, j->p, i -> Sin, 
raw_dat = squeeze(div_durs_compiled(:,:,j,:,:));

% i => Sin; j => p, 
lin_colors = lines(20);

fg = figure; hold on;
plot_indices = 1:length(secRatio);

xval = secRatio;
% fig 3e 3f 
for i = [2 7]
%	for j = 1:length(feed_rate)
	ct = 0;
	for j = plot_indices
		ct = ct + 1;
		doubDat = raw_dat(:,:,j,i)/60;
		doubDat = doubDat(:);
		doubDat = doubDat/mean(doubDat);
		
		dat_std(ct) = std(doubDat);
		dat_skew(ct) = skewness(doubDat);
		dat_kurt(ct) = kurtosis(doubDat);
		
		pd = fitdist(doubDat(:),'GeneralizedExtremeValue');
		gev_k(ct) = pd.k;
		gev_sig(ct) = pd.sigma;
		gev_mu(ct) = pd.mu;
		
	end

	subplot(6,1,1); hold on;
	plot(xval, dat_std,'-o','linewidth',2)
	xticks(xval)
	fg.Children(1).XScale = 'log';
	%ylim([0.43 0.47])
	ylabel('CV','FontSize',11)
	
	subplot(6,1,2); hold on;
	plot(xval,dat_skew,'-o','linewidth',2)
	ylabel('Skewness','FontSize',11)
	xticks(xval)
	fg.Children(1).XScale = 'log';
	
	subplot(6,1,3); hold on;
	plot(xval,dat_kurt,'-o','linewidth',2)
	ylabel('Kurtosis','FontSize',11)
	xticks(xval)
	fg.Children(1).XScale = 'log';
	
	subplot(6,1,4); hold on;
	plot(xval,gev_k,'-o','linewidth',2)
	ylabel('Shape (k)','FontSize',11)
	xticks(xval)
	fg.Children(1).XScale = 'log';
	
	subplot(6,1,5); hold on;
	plot(xval,gev_sig,'-o','linewidth',2)
	ylabel('Scale (σ)','FontSize',11)
	xticks(xval)
	fg.Children(1).XScale = 'log';
	
	subplot(6,1,6); hold on;
	plot(xval,gev_mu,'-o','linewidth',2)
	ylabel('Location (μ)','FontSize',11)
	xticks(xval)
	fg.Children(1).XScale = 'log';
	
	xlabel('Sec Ratio (Log)','FontSize',11)
end

fg.Position = [378	-11		275  	700];

name = 'fig3ef_CV_skew_kurt_logGEVpars_red_high_blue_low_Sin';


%figure2pdf(strcat(name,'.pdf'));
export_fig(name, '-png', '-r300', '-transparent', '-painters')


% ========================================

%%ANOVA	Fig 3g
%Sin_indices = 2:end;
% div_durs_compiled(:,:,j,l,i)		=> l - secRatio, j->p, i -> Sin, 
j = 1; 	% j = 1 => p = 3. Only p that data extracted
raw_dat = squeeze(div_durs_compiled(:,:,j,:,2:end));

dims = size(raw_dat);
%dims = dims([1 2 4 3]);
clear data
data = nan(prod(dims(1:3)),dims(4));

ct = 0;
for i = 1:dims(3)		
	for k = 1:dims(2)
		for l = 1:dims(1)
			ct = ct + 1;
			data(ct,:) = squeeze(raw_dat(l,k,:,i));
			%data(ct,:) = squeeze(raw_dat(l,k,i,:));
		end
	end
end


[p1,tbl,stats] = anova2(data,prod(dims(1:2)))

name = 'fig3g_anova';
writecell(tbl,strcat(name,'.csv'))

results = multcompare(stats,'Dimension',[1 2])
%f = figure(gcf)
%f.Position(3:4) = [460 135];
%f.Position(3:4) = [530 147];


%export_fig(name, '-png', '-r300', '-transparent', '-painters')