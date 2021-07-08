%% Figure 4 - Effect of Nutrient Uptake and Secretion

% Add Path for figure2pdf
addpath('C:\Users\sysbio_admin\Google Drive\MATLAB File Exchange\')

% ==============
%% (A) Prototrophs + Import of limiting nutrient. 
% ==============
p = 3;
seed_range = 0:2;
n_ext_range = [1 p];
feed_rate = [0 312 625 1250:1250:5000 7500 12500];

Sin = 1e3*[1.5 5 6 6.5 7 9 13];

% ==============
%% (i) Shifted Histogram, p=3, high Sin=50k,and low sin= 7.5k
% ==============
base_dir = 'D:\Debu Simulations\Sep 2020\var_extFeed\';

load(strcat(base_dir,'var_extFeed_fig4.mat'))

lin_colors = lines(14);
% div_durs_compiled(:,k,j,n_ext,i). i=> Sin, k => Seed, n_ext => No of metabs imported
n_ext = 1;

f = figure; 
%plot_indices = 1:2:9;
plot_indices = [1 2 3 4 9];
subplot(2,2,1); hold on;
%i = 1; 	% Sin = 1.5e3
i = 2; % Sin = 5e3
%for j = 1:length(feed_rate)
for j = 1:length(plot_indices) %1:length(feed_rate)
	c2(:,:) = div_durs_compiled(:,:,j,n_ext,i)/60;
	
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Feed =',num2str(feed_rate(plot_indices(j)))},' '));
	%h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Feed =',num2str(feed_rate(j))},' '));
end
%xlim([90 145])
xlim([25 55])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
%legend('location','northeast')

subplot(2,2,2); hold on;
%i = 5; 
i = 7; 	% Sin = 13e3
%for j = 1:length(feed_rate)
for j = 1:length(plot_indices) %1:length(feed_rate)
	c2(:,:) = div_durs_compiled(:,:,j,n_ext,i)/60;
	
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Feed =',num2str(feed_rate(plot_indices(j)))},' '));
	%h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Feed =',num2str(feed_rate(j))},' '));
end
xlim([20 45])
xlabel('Cell generation times (in mins)')
ylabel('PDF')

legend('location','northeast')	% ADJUST MANUALLY


%sgtitle({'Histograms of generation times, when one limiting metabolite is', 'fed directly at varying rates, and varying S_{in}'},'FontSize',10,'Fontweight','bold');

name = 'fig4a_1feedrate_distCompare_lowSin_highSin';
set(f, 'Color', 'none');
export_fig(name, '-png', '-r300', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));

% ==============
%% (ii) Growth rate vs Sin, for Different Metabolite Feed rates
% ==============
%% load original prototroph var_sin_p data
path_loc = 'C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Monod\';

ld_dat = load(strcat(path_loc,'divD_growR_vary_sin_p_sizer.mat'),'growth_rate');%,'Sin','p','seed_range');
growth_rate_o = ld_dat.growth_rate;
clear('ld_dat');
load(strcat(path_loc,'divD_growR_vary_sin_p_sizer.mat'),'x_dat');%,'Sin','p','seed_range');

%% Select the data points for manuscript Plot figures for Sin
sel = [4 6 7 8 10 13];
% growth_rate_o(k,i,j) => k-Seed, i - Sin, j - feed_rate

% growth_rate(k,j,n_ext,i) 
plot_indices = [1 2 3 4 9];
n_ext = 1;
fg = figure; hold on;
%% Plot the p = 2 dataset. j = 2 => p=2
%plot(Sin, reshape(mean(growth_rate_o(:,:,2),1),length(Sin),[]),'--o','Displayname','p = 2','linewidth',2,'color','k')
plot(Sin(2:end), reshape(mean(growth_rate_o(:,sel,2),1),length(Sin(2:end)),[]),'--o','Displayname','p = 2','linewidth',2,'color','k')

% plot p = 3 equivalent to feed = 0; with lines
j = 1;
plot(Sin(2:end), reshape(mean(growth_rate(:,j,n_ext,2:end),1),length(Sin(2:end)),[]),'-o','Displayname',strjoin({'Feed =',num2str(feed_rate(plot_indices(j))),'(p = 3)'},' '),'linewidth',2,'color',lin_colors(j,:))

for j = 2:length(feed_rate(plot_indices))
	plot(Sin(2:end), reshape(mean(growth_rate(:,j,n_ext,2:end),1),length(Sin(2:end)),[]),'-o','Displayname',strjoin({'Feed =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',2,'color',lin_colors(j,:))
end
legend('location','southeast','NumColumns',2);
%title({'Variation of growth rate with Substrate flux','and limiting metabolite feed rate',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux available per cascade (molecules per sec)')
xlim([4800 13100])

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = 'fig4a_growthR_feedR';
set(fg, 'Color', 'none');
export_fig(name, '-png', '-pdf','-r300', '-transparent', '-painters')
%figure2pdf(strcat(name,'.pdf'));

%savefig(strcat(name,'.fig'));
%print(name,'-dpng','-r300');

%=====================
% div_durs_compiled(:,k,j,n_ext,i) 
%c2(:,:) = div_durs_compiled(:,:,j,n_ext,i)/60;

n_ext = 1;
figure; hold on;
%% Plot the p = 2 dataset. j = 2 => p=2
%plot(Sin, reshape(mean(growth_rate_o(:,:,2),1),length(Sin),[]),'--o','Displayname','p = 2','linewidth',2,'color','k')
%plot(Sin(2:end), reshape(mean(growth_rate_o(:,sel,2),1),length(Sin(2:end)),[]),'--o','Displayname','p = 2','linewidth',2,'color','k')
for j = 1:length(feed_rate)
	c2 = permute(div_durs_compiled(:,:,j,n_ext,:)/60,[5 1 2 3 4]);
	plot(Sin', mean(mean(c2,2),3), '-o','Displayname',strjoin({'Feed =',num2str(feed_rate(j))},' '),'linewidth',2)
end
legend('location','southeast','NumColumns',2);
title({'Variation of growth rate with Substrate flux','and limiting metabolite feed rate',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux available per cascade (molecules per sec)')

name = 'fig4a_growthR_feedR';
savefig(strcat(name,'.fig'));
print(name,'-dpng','-r300');





% ==============
%% (iii) Metabolite Uptake vs Change in growth rate
% ==============

plot_indices = [1:8];

% growth_rate(k,j,n_ext,i) 
fg = figure; hold on;
n_ext = 1;
for i = 2:length(Sin)
	%dat = (reshape(mean(growth_rate(:,plot_indices,n_ext,i),1),length(feed_rate(plot_indices)),[])-mean(growth_rate(:,1,n_ext,i),1))/mean(growth_rate(:,1,n_ext,i),1)*100;
	dat = ( mean(growth_rate(:,plot_indices,n_ext,i))./mean(growth_rate(:,1,n_ext,i)) - 1 )*100 ;
	% Only 3 seed values
%	dat = ( mean(growth_rate(1:3,plot_indices,n_ext,i))./mean(growth_rate(1:3,1,n_ext,i)) - 1 )*100 ;	
	plot(feed_rate(plot_indices), dat,'-o','Displayname',strjoin({'S_{in} =',num2str(Sin(i))},' '),'linewidth',2)
end
legend('location','southeast','NumColumns',2);
%title({'Change in growth rate with feed rate of limting metabolite for varying Substrate flux',''});
ylabel('% Change in growth rate')
xlabel('Metabolite feed rate (molecules per sec)')

%fg.Children(2).XScale = 'log'; %'linear'
fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = 'fig4a_PercentChange_growthR_feedR';
set(fg, 'Color', 'none');
export_fig(name, '-png', '-pdf','-r300', '-transparent', '-painters')


% Figure for REVIEWERS Comment regarding the lack of smoothness

%name = 'fig4a_PercentChange_growthR_feedR_reviewers_smooth';
name = 'fig4a_PercentChange_growthR_feedR_reviewers_rough';
%set(fg, 'Color', 'none');
export_fig(name, '-png', '-pdf','-r300', '-transparent', '-painters')


%savefig(strcat(name,'.fig'));
%print(name,'-dpng','-r300');


% =====================================================
%% (b) Prototroph + Secretion of nutrients.
% ==============
%% load data
base_dir = 'D:\Debu Simulations\Sep 2020\var_protSec\';
load(strcat(base_dir,'var_protSec_fig4.mat'),'div_durs_compiled','growth_rate');

Sin = 1e3*[1.5 5 6 6.5 7 9 13];
% p = [1 3 5];
p = 3;
seed_range = 0:2;
secRatio = [0 0.05 0.1 0.25 0.5 0.75];

% ==============
%% (i) Shifted Histogram. p=3, high Sin=50k, and low Sin = 7.5k
% ==============

lin_colors = lines(14);
% div_durs_compiled(:,k,j,l,i). i=> Sin, k => Seed, l => SecRatio

f = figure; 
subplot(2,2,1); hold on;
i = 2; 	% Sin = 5e3
j = 1; 	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	c2 = div_durs_compiled(:,:,j,l,i)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(l) = fill(h2,h1,lin_colors(l,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Sec. ratio =',num2str(secRatio(l))},' '));
end
xlim([25 170])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
%legend('location','northeast')

subplot(2,2,2); hold on;
i = 7; 
j = 1; 	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	c2 = div_durs_compiled(:,:,j,l,i)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(l) = fill(h2,h1,lin_colors(l,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Sec. ratio =',num2str(secRatio(l))},' '));
end
xlim([23 79])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
legend('location','northeast')

%sgtitle({'Histograms of generation times, when one limiting metabolite is', 'secreted at varying ratios of production, at varying S_{in}'},'FontSize',10,'Fontweight','bold')

name = 'fig4b_protSec_distCompare_lowSin_highSin';
set(f, 'Color', 'none');
export_fig(name, '-png', '-r300', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));

%savefig(strcat(name,'.fig'));
%print(name,'-dpng','-r300');

% ==============
% Standardized Distributions
% ==============

f = figure; 
hold on;
i = 2; 	% Sin = 5e3
j = 1; 	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	c2 = div_durs_compiled(:,:,j,l,i)/60;
%	c2 = (c2(:) - mean(c2(:)))/std(c2(:));
	c2 = c2(:)/mean(c2(:));
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(l) = fill(h2,h1,lin_colors(l,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Sec. ratio =',num2str(secRatio(l))},' '));
end
xlim([0.7 1.5])
ylim([0 6.5])
xlabel('Rescaled generation times')
ylabel('PDF')

f.Position(3:4) = f.Position(3:4)*0.6;
f.Children.XLabel.FontSize = 11;
f.Children.YLabel.FontSize = 11;

name = 'fig4b_protSec_STDdistCompare_lowSin';
set(f, 'Color', 'none');
export_fig(name, '-png', '-r300', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));
%print(name,'-dpng','-r300');


f = figure; hold on;
i = 7; 
j = 1; 	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	c2 = div_durs_compiled(:,:,j,l,i)/60;
%	c2 = (c2(:) - mean(c2(:)))/std(c2(:));
	c2 = c2(:)/mean(c2(:));
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(l) = fill(h2,h1,lin_colors(l,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Sec. ratio =',num2str(secRatio(l))},' '));
end
xlim([0.7 1.5])
ylim([0 6])
xlabel('Rescaled generation times')
ylabel('PDF')
%legend('location','northeast')
%title({'high Sin = 50k, Standardized generation time distribution, vary Secretion',''})
f.Position(3:4) = f.Position(3:4)*0.6;
f.Children.XLabel.FontSize = 11;
f.Children.YLabel.FontSize = 11;

name = 'fig4b_protSec_STDdistCompare_highSin';
set(f, 'Color', 'none');
export_fig(name, '-png', '-r300', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));
%print(name,'-dpng','-r300');

% ==============
%% (b) ii) Compare Growth rate, vs all Sin, for varying P_out.
% ==============

% growth_rate(k,j,l,i) 
f = figure; hold on;
j = 1; %=> p = 3;
for l = 1:length(secRatio)
	plot(Sin(2:end), reshape(mean(growth_rate(:,j,l,(2:end)),1),length(Sin)-1,[]),'-o','Displayname',strjoin({'Sec. ratio =',num2str(secRatio(l))},' '),'linewidth',2)
	%plot(Sin, reshape(mean(growth_rate(:,j,l,:),1),length(Sin),[]),'-o','Displayname',strjoin({'Sec. ratio =',num2str(secRatio(l))},' '),'linewidth',2)
	
%	c2 = permute(div_durs_compiled(:,:,j,l,:)/60,[5 1 2 3 4]);
%	plot(Sin', mean(mean(c2,2),3), '-o','Displayname',strjoin({'Sec. ratio =',num2str(secRatio(l))},' '),'linewidth',2)
	
	
end
legend('location','southeast')%,'NumColumns',3);
%title({'Variation of growth rate with Secretion ratio','varying Substrate flux',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux available per cascade (molecules per sec)')

f.Position(3:4) = f.Position(3:4)*0.75;
f.Children(2).XLabel.FontSize = 11;
f.Children(2).YLabel.FontSize = 11;
xl = xlim;
xlim([xl(1)-250 xl(2)+250])

name = 'fig4b_growthR_Sin_secRatio';
set(f, 'Color', 'none');
export_fig(name, '-png', '-pdf','-r300', '-transparent', '-painters')
%figure2pdf(strcat(name,'.pdf'));


% ==============
%% iii) Secretion vs Change in growth rate
% ==============

% growth_rate(k,j,l,i) 
f = figure; hold on;
j = 1; %=> p = 3;
for i = 1:length(Sin)
	dat = (reshape(mean(growth_rate(:,j,:,i),1),length(secRatio),[])-mean(growth_rate(:,j,1,i),1))/mean(growth_rate(:,j,1,i),1)*100;
	plot(secRatio, dat,'-o','Displayname',strjoin({'S_{in} =',num2str(Sin(i))},' '),'linewidth',2)
end
legend('location','southwest')%,'NumColumns',3);
%title({'Change in growth rate with Secretion ratio varying Substrate flux',''});
ylabel('% Change in growth rate')
xlabel('Secretion ratio')

f.Position(3:4) = f.Position(3:4)*0.75;
f.Children(2).XLabel.FontSize = 11;
f.Children(2).YLabel.FontSize = 11;

name = 'fig4b_PercentChange_growthR_secRatio';
set(f, 'Color', 'none');
export_fig(name, '-png', '-pdf','-r300', '-transparent', '-painters')




% ==============
%% Secration Ratio to Mean Secretion Rate
% ==============
2e7 * 0.1 = 2e6;
time = 50 mins = 3000 sec;
=> 2e6/3e3 = 0.67e3 molecules/sec = ~670/sec










%% =============================================================================================
%% Code for figures when No Secretion simulation not included in the simulation dataset

%% (b) Prototroph + Secretion of nutrients.
% ==============
%% load data
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';
load(strcat(base_dir,'\protSec\Synced\Sizer\','var_protSec_Dat_p3.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate');

Sin = 1e3*[5 7.5 11 20 50];
% p = [1 3 5];
p = 3;
seed_range = 0:2;
secRatio = [0.05 0.1 0.25 0.5 0.75];

% Load Dataset with no secretion
noSec = load('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\manuscript\Fig2\sin_p_var_sizerDat_new.mat')

% ==============
%% (i) Shifted Histogram. p=3, high Sin=50k, and low Sin = 7.5k
% ==============

lin_colors = lines(14);
% div_durs_compiled(:,k,j,l,i). i=> Sin, k => Seed, l => SecRatio
n_ext = 1;

f = figure; 
subplot(2,1,1); hold on;
% Plot noSec histogram
	i = 2; j = 3;	% Sin = 7.5k, p = 3
	c2 = noSec.x_dat(:,:,i,j)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(1) = fill(h2,h1,lin_colors(1,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName','No Secretion');
	
i = 1; j = 1; 	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	c2 = div_durs_compiled(:,:,j,l,i)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(l+1) = fill(h2,h1,lin_colors(l+1,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Sec. ratio =',num2str(secRatio(l),'%.2f')},' '));
end
xlim([0 600])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
%legend('location','northeast')

subplot(2,1,2); hold on;
% Plot noSec histogram
	i = 5; j = 3;	% Sin = 50k, p = 3
	c2 = noSec.x_dat(:,:,i,j)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(1) = fill(h2,h1,lin_colors(1,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName','No Secretion');
	
i = 5; j = 1; 	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	c2 = div_durs_compiled(:,:,j,l,i)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(l+1) = fill(h2,h1,lin_colors(l+1,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Sec. ratio =',num2str(secRatio(l),'%.2f')},' '));
end
xlim([0 250])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
legend('location','northeast')

sgtitle({'Histograms of generation times, when one limiting metabolite is', 'secreted at varying ratios of production, at varying S_{in}',''},'FontSize',10,'Fontweight','bold')

name = 'fig4b_protSec_distCompare_lowSin_highSin';
figure2pdf(strcat(name,'.pdf'));
print(name,'-dpng','-r300');

% ==============
% Standardized Distributions
% ==============
f = figure; 
hold on;
% Plot noSec histogram
	i = 2; j = 3;	% Sin = 7.5k, p = 3
	c2 = noSec.x_dat(:,:,i,j)/60;
	c2 = (c2(:) - mean(c2(:)))/std(c2(:));
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(1) = fill(h2,h1,lin_colors(1,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName','No Secretion');
	
i = 1; j = 1;	% j = 1 => p = 3. Only p that data extracted
for l = 1:length(secRatio)
	c2 = div_durs_compiled(:,:,j,l,i)/60;
	c2 = (c2(:) - mean(c2(:)))/std(c2(:));
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(l+1) = fill(h2,h1,lin_colors(l+1,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Sec. ratio =',num2str(secRatio(l),'%.2f')},' '));
end
xlim([-3 6])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
legend('location','northeast')
title({'Low Sin = 7.5k, Standardized generation time distribution, vary Secretion',''})

name = 'fig4b_protSec_STDdistCompare_lowSin';
figure2pdf(strcat(name,'.pdf'));
print(name,'-dpng','-r300');


figure; hold on;
% Plot noSec histogram
	i = 5; j = 3;	% Sin = 50k, p = 3
	c2 = noSec.x_dat(:,:,i,j)/60;
	c2 = (c2(:) - mean(c2(:)))/std(c2(:));
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(1) = fill(h2,h1,lin_colors(1,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName','No Secretion');
	
i = 5; j = 1;
for l = 1:length(secRatio)
	c2 = div_durs_compiled(:,:,j,l,i)/60;
	c2 = (c2(:) - mean(c2(:)))/std(c2(:));
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(l+1) = fill(h2,h1,lin_colors(l+1,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Sec. ratio =',num2str(secRatio(l),'%.2f')},' '));
end
xlim([-3 6])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
legend('location','northeast')
legend('location','northeast')
title({'high Sin = 50k, Standardized generation time distribution, vary Secretion',''})

name = 'fig4b_protSec_STDdistCompare_highSin';
figure2pdf(strcat(name,'.pdf'));
print(name,'-dpng','-r300');

% ==============
%% (b) ii) Compare Growth rate, vs all Sin, for varying P_out.
% ==============
% load noSec Data
noSec = load('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\manuscript\Fig2\growR_var_sin_p.mat','growth_rate');

% growth_rate(k,j,l,i) 
figure; hold on;
j = 3;
plot(Sin(2:end), reshape(mean(noSec.growth_rate(:,(2:end),j),1),length(Sin)-1,[]),'-o','Displayname','No Secretion','linewidth',2)

j = 1; %=> p = 3;
for l = 1:length(secRatio)
	plot(Sin(2:end), reshape(mean(growth_rate(:,j,l,(2:end)),1),length(Sin)-1,[]),'-o','Displayname',strjoin({'Sec. ratio =',num2str(secRatio(l),'%.2f')},' '),'linewidth',2)
end
legend('location','southeast')%,'NumColumns',3);
title({'Variation of growth rate with Secretion ratio','varying Substrate flux',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux available per cascade (molecules per sec)')


name = 'fig4b_growthR_Sin_secRatio';
figure2pdf(strcat(name,'.pdf'));
print(name,'-dpng','-r300');
