% Load data - SUPP for Fig 3a

p = 3;
seed_range = 0:6;
n_ext_range = [1 p];
feed_rate = [0 312 625 1250:1250:5000 7500 12500];

Sin = 1e3*[1.5 5 6 6.5 7 9 13];

base_dir = 'D:\Debu Simulations\Sep 2020\var_extFeed\';

load(strcat(base_dir,'var_extFeed_fig4_seed_0-6.mat'),'div_durs_compiled');

plot_indices = [1 2 3 4 7 9];

% ========================================
% STsysbio
ANDARDIZED => RESCALED TO MEAN

n_ext = 1;
%div_durs_compiled(:,k,j,n_ext,i)
raw_dat = squeeze(div_durs_compiled(:,:,:,n_ext,:));

% i => Sin; j => p, 
i = 5; j = 3;
% COLORS
lin_colors = lines(20);

fg = figure; hold on;

%for i = 1:length(Sin)
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
%	plot(feed_rate(plot_indices),dat_std,'-o','linewidth',2)
	plot(feed_rate(plot_indices),dat_std,'-o','linewidth',2)
	xticks(feed_rate(plot_indices))
	fg.Children(1).XScale = 'log';
	%ylim([0.43 0.47])
	ylabel('CV','FontSize',11)
	
	subplot(6,1,2); hold on;
	plot(feed_rate(plot_indices),dat_skew,'-o','linewidth',2)
	ylabel('Skewness','FontSize',11)
	xticks(feed_rate(plot_indices))
	fg.Children(1).XScale = 'log';
	
	subplot(6,1,3); hold on;
	plot(feed_rate(plot_indices),dat_kurt,'-o','linewidth',2)
	ylabel('Kurtosis','FontSize',11)
	xticks(feed_rate(plot_indices))
	fg.Children(1).XScale = 'log';
	
	subplot(6,1,4); hold on;
	plot(feed_rate(plot_indices),gev_k,'-o','linewidth',2)
	ylabel('Shape (k)','FontSize',11)
	xticks(feed_rate(plot_indices))
	fg.Children(1).XScale = 'log';
	
	subplot(6,1,5); hold on;
	plot(feed_rate(plot_indices),gev_sig,'-o','linewidth',2)
	ylabel('Scale (σ)','FontSize',11)
	xticks(feed_rate(plot_indices))
	fg.Children(1).XScale = 'log';
	
	subplot(6,1,6); hold on;
	plot(feed_rate(plot_indices),gev_mu,'-o','linewidth',2)
	ylabel('Location (μ)','FontSize',11)
	xticks(feed_rate(plot_indices))
	fg.Children(1).XScale = 'log';
	
	xlabel('Feed Rate (Log)','FontSize',11)
end

fg.Position = [378	-11		275  	700];

name = 'fig3ab_CV_skew_kurt_logGEVpars_red_high_blue_low_Sin';


%figure2pdf(strcat(name,'.pdf'));
export_fig(name, '-png', '-r300', '-transparent', '-painters')


% ========================================

%%ANOVA	Fig 3c
feed_indices = [1 2 3 4 9];
%Sin_indices = 2:end;
raw_dat = squeeze(div_durs_compiled(:,:,feed_indices,n_ext,2:end));

dims = size(raw_dat);
dims = dims([1 2 4 3]);
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

name = 'fig3c_anova';
writecell(tbl,strcat(name,'.txt'))
%f = figure(gcf)
%f.Position(3:4) = [460 135];
%f.Position(3:4) = [530 147];


%export_fig(name, '-png', '-r300', '-transparent', '-painters')


% ===============================================

% Overlap of p = 2 and prototroph feed = 7500

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
	
% Plot prototroph Sin = 13e3 (i = 7), feed = 7500 (l = 8)
n_ext = 1;
raw_dat = squeeze(div_durs_compiled(:,:,:,n_ext,:));		% div_durs_compiled(:,k,j,n_ext,i)

	j = 2; 	% colour and order of graph
	l = 8; i = 7;
	c2 = raw_dat(:,:,l,i)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Prot. Fd =',num2str(feed_rate(l)),'CV =',num2str(std(c2(:)/mean(c2(:))),'%.2g')},' '));

xlabel('Cell generation times (in mins)')
ylabel('PDF')
legend('location','northeast','NumColumns',1)

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = strcat('Hist_compare_p2_sin13k_ExtFeed7500');
set(fg, 'Color', 'none');
export_fig(name, '-png', '-r300', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));

