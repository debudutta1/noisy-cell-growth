%% Figure 6 - Cross-feeding

%% LOAD DATA

%% CODE WORKS ON MY LAB PC and not LAB WORKSTATION <<=============

% Load data
sec_or_min = 0; % sec
n = 3;
p = 3;
n_ext = 0;
seed_range = 0:3;
Sin = 1e3*[1.5 5 6 6.5 7 9 13];
feed_rate = [312 625 1250 2500 5000 7500];% 12500];
secRatio = [0.3:0.1:0.7];
gens = 13; % 2^14 cells = 16384
base_dir = 'D:\Debu Simulations\Sep 2020\var_ALT_auxFeedSec\1_fluxDouble\';
load(strcat(base_dir,'1_fluxDouble_ALT_Dat.mat'),'div_durs_compiled','sec_end_compiled','growth_rate');

%=============== VERSION 2 of DATA
% Load data
sec_or_min = 0; % sec
n = 3;
p = 3;
n_ext = 0;
seed_range = 0:3;
Sin = 1e3*[4 5 6 6.5 7 9];
feed_rate = 1e3*[4 5 5.5 6 6.5 7];% 12500];
secRatio = [0.3:0.1:0.7];

base_dir = 'D:\Debu Simulations\Sep 2020\var_ALT_auxFeedSec\1_fluxDouble\';
load(strcat(base_dir,'1_fluxDouble_ALT_Dat_v2.mat'),'div_durs_compiled','sec_end_compiled','growth_rate');

%============

%nofeed = load('C:\Users\Dibyendu_Lab\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\manuscript\Fig2\growR_var_sin_p.mat');
%nofeed.Sin = [5000 7500 11000 20000 50000];

nofeed = load('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Monod\divD_growR_vary_sin_p_Sizer.mat');

% ExtFeed Prototroph data with Feed
extfeed = load(strcat('D:\Debu Simulations\Sep 2020\var_extFeed\var_extFeed_fig4.mat'));

% AuxFeed
auxfeed = load('D:\Debu Simulations\Sep 2020\var_auxFeed\auxFeed_fig5_dat2.mat','div_durs_compiled','growth_rate'); 
auxfeed.Sin = 1e3*[1.5 3 4 5 6 6.5 7 9 13];

% =========================
% NEW FIGURE - 5b: Growth rate vs Sin. No Feed, ExtFeed, AuxFeed, AuxFeedSec
% =========================

%sel = [2 3 4 6 7 8 10 13];	% %% Select the data points for Sin - nofeed

% COLORS
lin_colors = lines(20);

% growth_rate(k,j,i) => k-Seed, j - feed_rate, i - Sin
f = figure; hold on;

% Plot noFeed
sel = [3 4 6 7 8 10];	% %% Select the data points for Sin - nofeed
j = 3; %=>p=3;
plot(nofeed.Sin(sel), reshape(mean(nofeed.growth_rate(:,sel,j),1),length(nofeed.Sin(sel)),[]),'--o','Displayname','Prototroph','linewidth',2,'color',[.4 .4 .4])

% Plot auxFeed
j = 6; %=>feed_rate = 7500 => Saturated effect;
sel2 = [3:8];
plot(auxfeed.Sin(sel2), reshape(mean(auxfeed.growth_rate(:,j,sel2),1),length(auxfeed.Sin(sel2)),[]),'--o','Displayname','Auxotroph','linewidth',2,'color','k')


sel_dat = 1:length(Sin);
%plot_indices = 4:7;	% feed_rate
plot_indices = 1:5;	% feed_rate
l = 3;% =>0.5
%for j = 2:length(feed_rate)
for j = 1:length(plot_indices)
	plot(Sin(sel_dat), reshape(mean(growth_rate(:,plot_indices(j),l,sel_dat),1),length(Sin(sel_dat)),[]),'-o','Displayname',strjoin({'Feed =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',2,'color',lin_colors(j,:))
end
%title({'Variation of growth rate of auxotroph with varying feed','and varying Substrate flux',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux (per sec)')
%legend('location','east')%,'NumColumns',3);
%legend('location','southoutside','NumColumns',3);
legend('location','east','NumColumns',2);

f.Position(3:4) = f.Position(3:4)*0.7;
f.Children(2).XLabel.FontSize = 11;
f.Children(2).YLabel.FontSize = 11;

name = 'fig6_new_noFeed_auxFeed_auxFeedSec';
set(f, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r300', '-transparent', '-painters')


for j = plot_indices
	% Plot extFeed => % growth_rate = nan(length(seed_range),length(feed_rate),length(n_ext_range),length(Sin));
	plot(Sin(sel_dat), reshape(mean(extFeed.growth_rate(:,plot_indices(j),1,sel),1),length(Sin(sel_dat)),[]),'--o','Displayname',strjoin({'p =', num2str(j)},' '),'linewidth',2,'color','k')
	
	% PLOT AuxFeedSec => %% growth_rate(k,j,l,i)=> k-Seed, j-feed_rate, l-secRatio, i-Sin, 
	plot(Sin(sel_dat), reshape(mean(growth_rate(:,plot_indices(j),l,sel_dat),1),length(Sin(sel_dat)),[]),'-o','Displayname',strjoin({'Feed =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',2,'color',lin_colors(j,:))
	
end






% =========================

%% Fig -6 (c) : Bubble colouring based on growth rate faster than prototroph no feed.

% Plot X: Feed_rate, Y: Mean Secretion Rate. Scatter with bubbles of sizes based on growth rate. Fill bubbles that grow faster than the prototroph
lin_colors = lines(14);
%% (:,k,j,l,i)=> 1st-all cells, k-Seed, j-feed_rate, l-secRatio, i-Sin, 

% High Sin
i = 4; 
f = figure; hold on;

plot_indices = [1 3 6];	%for Sin
sel = [3 6 10];
j = 3;
%for i = 1:length(Sin)
for i = 1:length(plot_indices)
	f = figure; hold on;

	cdat = permute(nanmean(growth_rate(:,:,:,plot_indices(i))),[3 2 1 4]);
	for l = 1:length(secRatio)
		% Mean Secretion rate from Secretion Ratio
		y = permute(mean(mean(sec_end_compiled(:,:,:,l,plot_indices(i))./div_durs_compiled(:,:,:,l,plot_indices(i)))),[4 3 2 1]);
		
		% Marker size based on growth_rate
		sz = fix((cdat(l,:) - nanmin(cdat(:)))/(nanmax(cdat(:)) - nanmin(cdat(:)))*100)+35;
		sz = fix((cdat(l,:) - nanmin(cdat(:)))/(nanmax(cdat(:)) - nanmin(cdat(:)))*100)+50;
		
		h(l,i) = scatter(feed_rate, y,sz,'MarkerEdgeColor',lin_colors(i,:));
		
		%select = cdat(l,:) > round(max(mean(nofeed.growth_rate(:,:,1,i))),2);
		select = cdat(l,:) > mean(nofeed.growth_rate(:,sel(i),j));
		if sum(select) > 0
			fl(l,i) = scatter(feed_rate(select), y(select), sz(select),'MarkerFaceColor',lin_colors(i,:),'MarkerEdgeColor',lin_colors(i,:),'LineWidth',0.7);
		end
	end
	leg = scatter(nan,nan,'MarkerFaceColor',lin_colors(i,:),'MarkerEdgeColor',lin_colors(i,:),'Displayname',strjoin({'S_{in} =',num2str(Sin(plot_indices(i)))},' '));
	%leg(i) = scatter(nan,nan,'MarkerFaceColor',lin_colors(i,:),'MarkerEdgeColor',lin_colors(i,:),'Displayname',strjoin({'S_{in} =',num2str(Sin(i))},' '));
	%pause
%end

% Make figure xlim ylim square
%xl = xlim;		yl = ylim;
%xlim([min(xl(1),yl(1)) max(xl(2),yl(2))])
%ylim([min(xl(1),yl(1)) max(xl(2),yl(2))])

% diagonal guideline
g = plot(xlim,xlim,'--','color',[.5 .5 .5],'linewidth',2);

%xticks(feed_rate)

%title({'Auxotrophy Advantage over prototroph',''});
ylabel('Mean Secretion rate (molecules per sec)')
xlabel('Metabolite feed rate (molecules per sec)')
legend(leg,'location','northeast');
%legend(leg,'location','southoutside','NumColumns',4);

%f.Position(2) = 100;
%f.Position(4) = 500;

f.Position(3:4) = f.Position(3:4)*0.7;
f.Children(2).XLabel.FontSize = 11;
f.Children(2).YLabel.FontSize = 11;
pause
name = strcat('fig6a_auxAdv_nofeed_sin',num2str(Sin(plot_indices(i))));
set(f, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r300', '-transparent', '-painters')

end



%savefig(strcat(name,'.fig'));
%print(name,'-dpng','-r600');
% =========================

