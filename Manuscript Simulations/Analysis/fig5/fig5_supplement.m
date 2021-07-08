% Load data - SUPP for Fig 5b

n = 3;
p = 3;
n_ext = 0;
seed_range = 0:3;
Sin = 1e3*[4 5 6 6.5 7 9];
feed_rate = 1e3*[4 5 5.5 6 6.5 7];% 12500];
secRatio = [0.3:0.1:0.7];

base_dir = 'D:\Debu Simulations\Sep 2020\var_ALT_auxFeedSec\1_fluxDouble\';
load(strcat(base_dir,'1_fluxDouble_ALT_Dat_v2.mat'),'div_durs_compiled');%,'sec_end_compiled','growth_rate');

lin_colors = lines(14);

% ========================================
% STANDARDIZED => RESCALED TO MEAN

%i = 2	%low sin
%i = 7	% high sin

clear raw_dat doubDat doubDat dat_skew dat_kurt gev_k gev_sig gev_mu dat_std

% i => Sin; j => p, 
lin_colors = lines(20);

fg = figure; hold on;

plot_indices = 4:7;	% feed_rate
%plot_indices = 4:6;		% Since growth rates of 6,7 are identical in effect; feed_rate
xval = feed_rate(plot_indices);

% div_durs_compiled(:,:,l,i)		=> l - feed_rate, i -> Sin, 
raw_dat = squeeze(div_durs_compiled(:,:,plot_indices,:));

% fig 4b
for i = [2 7]
%	for j = 1:length(feed_rate)
	ct = 0;
	for j = 1:length(plot_indices)
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
	
	xlabel('Feed Rate (Log)','FontSize',11)
end

fg.Position = [378	-11		275  	700];

name = 'fig4b_CV_skew_kurt_logGEVpars_red_high13k_blue_low3k_Sin';


%figure2pdf(strcat(name,'.pdf'));
export_fig(name, '-png', '-r300', '-transparent', '-painters')


% ========================================

%%ANOVA	Fig 5b
%Sin_indices = 2:end;
% div_durs_compiled(:,:,l,i)		=> l - feed_rate, i -> Sin, 
clear data raw_dat

plot_indices = 1:5;		% feed_rate

raw_dat = squeeze(div_durs_compiled(:,:,plot_indices,2:end));

dims = size(raw_dat);
%dims = dims([1 2 4 3]);
data = nan(prod(dims(1:3)),dims(4));

ct = 0;
for i = 1:dims(3)		
	for k = 1:dims(2)
		for l = 1:dims(1)
			ct = ct + 1;
			%data(ct,:) = squeeze(raw_dat(l,k,:,i));
			data(ct,:) = squeeze(raw_dat(l,k,i,:));
		end
	end
end


[p1,tbl,stats] = anova2(data,prod(dims(1:2)))

name = 'fig5b_anova';
writecell(tbl,strcat(name,'.csv'))

results = multcompare(stats,'Dimension',[1 2])
%f = figure(gcf)
%f.Position(3:4) = [460 135];
%f.Position(3:4) = [530 147];


%export_fig(name, '-png', '-r300', '-transparent', '-painters')

% ===============================================

% Overlap of p = 2 and auxotroph feed = 7500

%% load var_sin_p data
Sin = 1e3*[13];
p = [2];
seed_range = 0:2;
% Compile Data from Raw Simulation data
base_dir = 'D:\Debu Simulations\Sep 2020\var_p_sin\';
% SIZER
fold_name = 'Sizer';
for j = 1:length(p)
	for i = 1:length(Sin)
		for k = seed_range
			clear div_durs div_durs_exp sim_vars
			load(strcat(base_dir,fold_name,'\var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs_exp');
			x_dat(k+1,:,i,j) = div_durs_exp;
		end
	end
end

fg = figure; hold on;

	j = 1;
	c2 = x_dat(:)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'p = 2','CV =',num2str(std(c2(:)/mean(c2(:))),'%.2g')},' '));
	
% Plot auxotroph Sin = 13e3 (i = 7), feed = 7500 (l = 6)
	j = 2; l = 6; i = 7;
	c2 = div_durs_compiled(:,:,l,i)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Aux. Fd =',num2str(feed_rate(l)),'CV =',num2str(std(c2(:)/mean(c2(:))),'%.2g')},' '));

xlabel('Cell generation times (in mins)')
ylabel('PDF')
legend('location','northeast','NumColumns',1)

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = strcat('Hist_compare_p2_sin13k_AuxFeed7500');
set(fg, 'Color', 'none');
export_fig(name, '-png', '-r300', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));

