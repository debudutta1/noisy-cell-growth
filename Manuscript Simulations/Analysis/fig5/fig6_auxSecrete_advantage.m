%% Figure 6 - Auxotrophic cells may produce a surplus of metabolites, while exhibiting a competitive growth rate compared to the prototroph
%%======================
% (a) Production doubled by gene expression parameter change
%%	(i) Bubbles of different size based on fitness. For different S_in, as subplots.
%%======================

%% CODE WORKS ON MY LAB PC and not LAB WORKSTATION <<=============

% Load data
sec_or_min = 0; % sec
n = 3;
p = 3;
n_ext = 0;
seed_range = 0:4;
Sin = 1e3*[7.5 10 13 20];
feed_rate = [7500:2500:15000];
secRatio = [0.3:0.1:0.7];
gens = 13; % 2^14 cells = 16384

addpath('C:\Users\Dibyendu_Lab\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\');

base_dir = 'E:\Debu_Work\Dec 2019\var_ALT_auxFeedSec\sizer\2_enzDouble_par7\';
load(strcat(base_dir,'2_enzDouble_par7_ALT_Dat.mat'),'div_durs_compiled','sec_end_compiled');
load(strcat(base_dir,'auxFeedSecrete_growR.mat'),'growth_rate');

base_dir = 'E:\Debu_Work\Dec 2019\var_ALT_auxFeedSec\sizer\1_fluxDouble\';
load(strcat(base_dir,'1_fluxDouble_ALT_Dat.mat'),'div_durs_compiled');
load(strcat(base_dir,'auxFeedSecrete_growR.mat'),'growth_rate');

%% (:,k,j,l,i)=> k-Seed, j-feed_rate, l-secRatio, i-Sin, 

% =========================
% Compute Mean Secretion Rate from Secretion Ratio

f = figure; hold on;
lin_colour = lines(length(Sin));
% Plot for each of the different cases mean secretion levels
for i = 1:length(Sin)
	for j = 1:length(feed_rate)
		h(i) = plot(secRatio, reshape(mean(mean(sec_end_compiled(:,:,j,:,i)./div_durs_compiled(:,:,j,:,i))),1,[]),'-o','linewidth', 2,'color',lin_colour(i,:),'displayname',strjoin({'S_{in}','=',num2str(Sin(i))},' '));
	end
end
xlabel('Secretion Ratio set')
ylabel('Resultant mean Secretion Rate')
title({'Secretion rate based on Secretion Ratio set',''})
legend(h,'location','northwest')
f.Children(2).YGrid = 'on';

% =========================
% ExtFeed Prototroph data with Feed

extfeed = load('C:\Users\Dibyendu_Lab\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\ExtFeed Prototroph New\var_extFeed_divdur_growR.mat');

%% growth_rate(k,j,n_ext,i)=> k-Seed, j-feed_rate, n_ext, i-Sin, (5,7,3,5)
% =========================
% Plot X: Feed_rate, Y: Mean Secretion Rate. Scatter with bubbles of sizes based on growth rate. Fill bubbles that grow faster than the prototroph
lin_colors = lines(14);
%% (:,k,j,l,i)=> 1st-all cells, k-Seed, j-feed_rate, l-secRatio, i-Sin, 

% High Sin
i = 4; 
ext_feed_sel = (4:7);
%ext_feed_sel = (4:-1:1);
f = figure; hold on;

for i = 1:length(Sin)
	cdat = permute(nanmean(growth_rate(:,:,:,i)),[3 2 1 4]);
	for l = 1:length(secRatio)
		% Mean Secretion rate from Secretion Ratio
		y = permute(mean(mean(sec_end_compiled(:,:,:,l,i)./div_durs_compiled(:,:,:,l,i))),[4 3 2 1]);
		
		% Marker size based on growth_rate
		sz = fix((cdat(l,:) - nanmin(cdat(:)))/(nanmax(cdat(:)) - nanmin(cdat(:)))*100)+35;
		
		h(l,i) = scatter(feed_rate, y,sz,'MarkerEdgeColor',lin_colors(i,:));
		
		%select = cdat(l,:) > round(max(mean(extfeed.growth_rate(:,:,1,i))),2);
		select = cdat(l,:) > mean(extfeed.growth_rate(:,ext_feed_sel,1,i));
		if sum(select) > 0
			fl(l,i) = scatter(feed_rate(select), y(select), sz(select),'MarkerFaceColor',lin_colors(i,:),'MarkerEdgeColor',lin_colors(i,:),'LineWidth',0.7);
		end
	end
	leg(i) = scatter(nan,nan,'MarkerFaceColor',lin_colors(i,:),'MarkerEdgeColor',lin_colors(i,:),'Displayname',strjoin({'S_{in} =',num2str(Sin(i))},' '));
	%pause
end
% diagonal guideline
g = plot(xlim,xlim,'--','color',[.5 .5 .5],'linewidth',2);
xticks(feed_rate)

title({'Auxotrophy Advantage',''});
ylabel('Mean Secretion rate (molecules per sec)')
xlabel('Metabolite feed rate (molecules per sec)')
legend(leg,'location','southoutside','NumColumns',4);

f.Position(2) = 100;
f.Position(4) = 500;

name = 'fig6_aux_advantage';
savefig(strcat(name,'.fig'));
print(name,'-dpng','-r600');

%% =====================
lin_colors = lines(14);
%% (:,k,j,l,i)=> 1st-all cells, k-Seed, j-feed_rate, l-secRatio, i-Sin, 

% High Sin
i = 4; 
ext_feed_sel = (4:7);
%ext_feed_sel = (4:-1:1);
f = figure; hold on;

for i = 1:length(Sin)
	a(i) = subplot(2,2,i);
	hold on;
	cdat = permute(nanmean(growth_rate(:,:,:,i)),[3 2 1 4]);
	for l = 1:length(secRatio)
		% Mean Secretion rate from Secretion Ratio
		y = permute(mean(mean(sec_end_compiled(:,:,:,l,i)./div_durs_compiled(:,:,:,l,i))),[4 3 2 1]);
		
		% Marker size based on growth_rate
		sz = fix((cdat(l,:) - nanmin(cdat(:)))/(nanmax(cdat(:)) - nanmin(cdat(:)))*100)+35;
		
		h(l,i) = scatter(feed_rate, y,sz,'MarkerEdgeColor',lin_colors(i,:));
		
		%select = cdat(l,:) > round(max(mean(extfeed.growth_rate(:,:,1,i))),2);
		select = cdat(l,:) > mean(extfeed.growth_rate(:,ext_feed_sel,1,i));
		if sum(select) > 0
			fl(l,i) = scatter(feed_rate(select), y(select), sz(select),'MarkerFaceColor',lin_colors(i,:),'MarkerEdgeColor',lin_colors(i,:),'LineWidth',0.7);
		end
	end
	leg(i) = scatter(nan,nan,'MarkerFaceColor',lin_colors(i,:),'MarkerEdgeColor',lin_colors(i,:),'Displayname',strjoin({'S_{in} =',num2str(Sin(i))},' '));
	%pause

	ylabel({'Mean Secretion rate (molecules per sec)'})
	xlabel({'Metabolite feed rate','(molecules per sec)'})	
	
	xl = [6000 16000];
	if i == 1
		xl = [6000 12500];
	end
	
	g = plot(xl,xl,'--','color',[.5 .5 .5],'linewidth',2);
	xticks(feed_rate)
	%xlim(xl);
	
	%yl = [0.4 2]*1e4;
	%ylim(yl);
end
sgtitle({'Auxotrophy Advantage',''},'FontSize',10,'Fontweight','bold')
%legend(leg,'location','southoutside','NumColumns',4);

f.Position = [ 680   337   738   641]

name = 'fig6_aux_advantage_subplot3';
savefig(strcat(name,'.fig'));
print(name,'-dpng','-r600');
%% =====================




mean(extfeed.growth_rate(:,ib1,1,ib2(i)))


perute(mean(extfeed.growth_rate(:,ext_feed_sel,1,i)))



%%=============

[C1, ia1, ib1] = intersect(feed_rate,extfeed.feed_rate);
[C2, ia2, ib2] = intersect(Sin,extfeed.Sin);

round(permute(mean(growth_rate(:,:,:,ia2(i))),[2 3 1])./mean(extfeed.growth_rate(:,ib1,1,ib2(i)))',2)












%%======================
%% SUPPLEMENTARY
%%======================
% (b) Production doubled using enzyme parameter
%%======================
%%======================