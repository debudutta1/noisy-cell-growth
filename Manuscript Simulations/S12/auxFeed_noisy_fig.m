% Manuscript Figure 6 - Plot Codes - Effect of Auxotrophy
%% =========================================================
n = 3;
p = 3;
seed_range = 0:3;
Sin = 1e3*[1.5 3 4 5 6 6.5 7 9 13];
% Noisy Feed rates are fractions to multiply with k3.
feed_rate = [0.5 0.8 0.9 1 1.1 1.2];

base_dir = 'D:\Debu Simulations\Sep 2020\var_auxFeed\';

%% Load data
load(strcat(base_dir,'auxFeed_fig5_dat.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate'); 

%% load var_sin_p data
%% load original prototroph var_sin_p data
path_loc = 'C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Monod\';
ld_dat = load(strcat(path_loc,'divD_growR_vary_sin_p_sizer.mat'),'growth_rate');%,'Sin','p','seed_range');
growth_rate_o = ld_dat.growth_rate;
clear('ld_dat');
load(strcat(path_loc,'divD_growR_vary_sin_p_sizer.mat'),'x_dat');%,'Sin','p','seed_range');


%====================
% Line COLORS
lin_colors = lines(14);
% ==============
%% a) Histogram of gen times, with low P_in (higher gen time), and high P_in (lower gen time)
% ==============
% div_durs_compiled(:,k,j,i) => k-Seed, i - Sin, j - feed_rate

plot_indices = 1:6;	% feed_rate

f = figure; 
subplot(2,2,1); hold on;
hold on;
i = 1;	% Low Sin
%for j = 1:length(feed_rate)
for j = 1:length(plot_indices)
	c2 = div_durs_compiled(:,:,plot_indices(j),i)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Feed rate =',num2str(feed_rate(plot_indices(j)))},' '));
end
xlim([65 105])
ylim([0 0.11])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
%legend('location','northeast')

subplot(2,2,2); hold on;
i = 5;	% High Sin
%for j = 1:length(feed_rate)
for j = 1:length(plot_indices)
	c2 = div_durs_compiled(:,:,plot_indices(j),i)/60;
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(j) = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'Feed rate =',num2str(feed_rate(plot_indices(j)))},' '));
end
xlim([20 50])
ylim([0 0.2])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
legend('location','northeast','NumColumns',1)

sgtitle({'Histograms of generation times, when auxotroph is fed the', 'missing metabolite at varying S_{in}',''},'FontSize',10,'Fontweight','bold')

% Manually alter spacing inside the figure and axes to create space for legend
f.Children(3).Position(4) =f.Children(3).Position(4)*.8
f.Children(4).Position(4) =f.Children(4).Position(4)*.8

name = 'fig5a_auxfeed_distCompare_lowSin_highSin';
%name = 'fig5a_auxfeed_distCompare_lowSin_highSin_noLegend';
savefig(strcat(name,'.fig'));
print(name,'-dpng','-r300');


% ==============
%% (b) Growth rate vs Sin, prototroph vs auxotroph, low and high P_in
% ==============

%noSec = load('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\manuscript\Fig2\growR_var_sin_p.mat','growth_rate');
%plot(Sin, reshape(mean(noSec.growth_rate(:,(2:end),j),1),length(Sin),[]),'--o','Displayname','p = 3','linewidth',2,'color','k')

sel = [2 3 4 6 7 8 10 13];	% %% Select the data points for Sin - growth_rate_o
sel_dat = 2:length(Sin);

plot_indices = 1:6;	% feed_rate

% growth_rate(k,j,i) => k-Seed, j - feed_rate, i - Sin
fg = figure; hold on;
% Plot in black dotted growth rates of p = 2, and p = 3 prototroph	
j = 3; %=>p=3;
plot(Sin(sel_dat), reshape(mean(growth_rate_o(:,sel,j),1),length(Sin(sel_dat)),[]),'--o','Displayname',strjoin({'p =', num2str(j),'\newline'},' '),'linewidth',2,'color',[.4 .4 .4])
j = 2; %=>p=2;
plot(Sin(sel_dat), reshape(mean(growth_rate_o(:,sel,j),1),length(Sin(sel_dat)),[]),'--o','Displayname',strjoin({'p =', num2str(j)},' '),'linewidth',2,'color','k')

%for j = 2:length(feed_rate)
for j = 1:length(plot_indices)
	plot(Sin(sel_dat), reshape(mean(growth_rate(:,plot_indices(j),sel_dat),1),length(Sin(sel_dat)),[]),'-o','Displayname',strjoin({'Feed =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',2,'color',lin_colors(j,:))
end
%title({'Variation of growth rate of auxotroph with varying feed','and varying Substrate flux',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux (per sec)')
%legend('location','east')%,'NumColumns',3);
%legend('location','southoutside','NumColumns',3);
legend('location','east','NumColumns',1);

xlim([2700 13300])
ylim([0.57 1.37])
%xlim([1250 13250])

%addpath('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\MATLAB file exchange\Break Axis\breakyaxis\')
%breakyaxis([.33 .58]);

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = 'fig5b_growthR_Sin_feedR';
%savefig(strcat(name,'.fig'));
%print(name,'-dpng','-r300');
set(fg, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r300', '-transparent', '-painters')

% ==============
%% (c) Percentage change in growth due to auxotrophic feed
% ==============
plot_indices = 1:6;	% feed_rate
sel = [1 2 3 4 6 7 8 10 13];	% %% Select the data points for Sin - growth_rate_o
sel_dat = [2:6 9];%2:length(Sin);

% growth_rate(k,j,i) => k-Seed, j - feed_rate, i - Sin
fg = figure; hold on;
%plot([2e3 16e3],[0 0],'--','color','k','linewidth',2,'DisplayName','Prototroph')
plot([0 13.5e3],[0 0],'--','color','k','linewidth',2,'DisplayName','Prototroph')
cnt = 1;
%for i = 1:length(Sin)
for i = sel_dat
	j = 3; % => p = 3
%	dat = (reshape(mean(growth_rate(:,plot_indices,i),1),length(plot_indices),[]) - mean(growth_rate_o(:,:,j),1))/mean(growth_rate_o(:,:,j),1)*100;
	dat = (reshape(mean(growth_rate(:,plot_indices,i),1),length(plot_indices),[])/mean(growth_rate_o(:,sel(i),j),1) - 1)*100;
	plot(feed_rate(plot_indices), dat,'-o','Displayname',strjoin({'S_{in} =',num2str(Sin(i))},' '),'linewidth',2,'color',lin_colors(cnt,:))
	cnt = cnt + 1;
end
%title({'Variation of growth rate of auxotroph with varying feed','and varying Substrate flux',''});
ylabel('% Change in growth rate')
xlabel('Metabolite feed rate (molecules per sec)')
legend('location','southeast','NumColumns',2);

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = 'fig5c_PercentChange_growthR_feedR_Sin';
set(fg, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r300', '-transparent', '-painters')

% ==============
%% (d) Change in Growth rate reduces with increasing p.
% ==============
proto = load('fig5d_dat_sizer.mat','growth_rate','p','Sin');
%p = [1, 2, 3, 4, 5, 10, 11, 14, 15, 19, 20, 24, 25];
%Sin = 1e3*[3 6 7 13];

% Plot change from p => p-1, for [2 3 4 5 11 15 20 25]
% proto.growth_rate(k,i,j)
fg = figure; hold on;
selected = [2 3 4 5 7 9 11 13];	% selected p
%h = plot([0 25],[0 0],'--','color','k','linewidth',2);
%set( get( get( h, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );	% HIDE LEGEND ENTRY
for i = 1:length(proto.Sin)
	dat1 = reshape(mean(proto.growth_rate(:,i,selected),1),[],1);
	dat2 = reshape(mean(proto.growth_rate(:,i,selected-1),1),[],1);
	dat = (dat2./dat1 - 1)*100;
	plot(proto.p(selected), dat,'-o','Displayname',strjoin({'Sin =',num2str(proto.Sin(i))},' '),'linewidth',2)
	%bar(categorical(p(2:end)), -dat)%,'-o','Displayname',strjoin({'Sin =',num2str(Sin(i))},' '),'linewidth',2)
end
xticks(p(selected))
legend('location','northeast')%,'NumColumns',3);
%title({'% maximum increase in growth rate upon loss of a bottleneck pathway', 'for varying Substrate flux',''});
ylabel({'% increase in growth rate for','p → p - 1'})
xlabel('Number of bottleneck pathways (p)')

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = 'fig5d_PercentChange_growthR_var_p';
set(fg, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r300', '-transparent', '-painters')

% ==============
% ==============
% ==============