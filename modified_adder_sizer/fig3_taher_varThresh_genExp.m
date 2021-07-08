%% Fig 3 graphs

%% Taher data analysis

load('C:\Users\Dibyendu_Lab\Google Drive\Conference and Schools\CCCP 2020\MotherMachine data\Suckjoon Jun Lab\Taher2015_data.mat')

data{1}(:,1) = glycerol.gen_time_min_ss(glycerol.elongation_rate_ss > 0.01);
data{1}(:,2) = glycerol.s_d_um_ss(glycerol.elongation_rate_ss > 0.01);
data{1}(:,3) = glycerol.s_b_um_ss(glycerol.elongation_rate_ss > 0.01);
data{1}(:,4) = glycerol.elongation_rate_ss(glycerol.elongation_rate_ss > 0.01);

data{2}(:,1) = sorbitol.gen_time_min_ss(sorbitol.elongation_rate_ss > 0.01);
data{2}(:,2) = sorbitol.s_d_um_ss(sorbitol.elongation_rate_ss > 0.01);
data{2}(:,3) = sorbitol.s_b_um_ss(sorbitol.elongation_rate_ss > 0.01);
data{2}(:,4) = sorbitol.elongation_rate_ss(sorbitol.elongation_rate_ss > 0.01);

data{3}(:,1) = glucose.gen_time_min_ss(glucose.elongation_rate_ss > 0.01);
data{3}(:,2) = glucose.s_d_um_ss(glucose.elongation_rate_ss > 0.01);
data{3}(:,3) = glucose.s_b_um_ss(glucose.elongation_rate_ss > 0.01);
data{3}(:,4) = glucose.elongation_rate_ss(glucose.elongation_rate_ss > 0.01);

data{4}(:,1) = glucose6aa.gen_time_min_ss(glucose6aa.elongation_rate_ss > 0.01);
data{4}(:,2) = glucose6aa.s_d_um_ss(glucose6aa.elongation_rate_ss > 0.01);
data{4}(:,3) = glucose6aa.s_b_um_ss(glucose6aa.elongation_rate_ss > 0.01);
data{4}(:,4) = glucose6aa.elongation_rate_ss(glucose6aa.elongation_rate_ss > 0.01);

data{5}(:,1) = glucose12aa.gen_time_min_ss(glucose12aa.elongation_rate_ss > 0.01);
data{5}(:,2) = glucose12aa.s_d_um_ss(glucose12aa.elongation_rate_ss > 0.01);
data{5}(:,3) = glucose12aa.s_b_um_ss(glucose12aa.elongation_rate_ss > 0.01);
data{5}(:,4) = glucose12aa.elongation_rate_ss(glucose12aa.elongation_rate_ss > 0.01);

data{6}(:,1) = syntheticrich.gen_time_min_ss(syntheticrich.elongation_rate_ss > 0.01);
data{6}(:,2) = syntheticrich.s_d_um_ss(syntheticrich.elongation_rate_ss > 0.01);
data{6}(:,3) = syntheticrich.s_b_um_ss(syntheticrich.elongation_rate_ss > 0.01);
data{6}(:,4) = syntheticrich.elongation_rate_ss(syntheticrich.elongation_rate_ss > 0.01);

data{7}(:,1) = TSB.gen_time_min_ss(TSB.elongation_rate_ss > 0.01);
data{7}(:,2) = TSB.s_d_um_ss(TSB.elongation_rate_ss > 0.01);
data{7}(:,3) = TSB.s_b_um_ss(TSB.elongation_rate_ss > 0.01);
data{7}(:,4) = TSB.elongation_rate_ss(TSB.elongation_rate_ss > 0.01);

cond_list = {'Glycerol','Sorbitol','Glucose','Glucose+6AA','Glucose+12AA','Synthetic Rich','TSB'};


%%======================

% Additional filtering to remove zero generation time
for i = 1:7
	idx = find(data{i}(:,1)==0);
	data{i}(idx,:) = [];
	idx = find(data{i}(:,2)==data{i}(:,3));
	data{i}(idx,:) = [];
end

%%======================
% Plot histogram of inter-division times
%%======================
lin_colors = lines(7);

figure; hold on;
for i = 1:7
	subplot(7,1,i);
	hold on;
	if i < 3
		histogram(data{i}(:,1),'normalization','pdf','BinMethod','fd','BinWidth',4,'EdgeColor','none','Displayname',cond_list{i},'FaceColor',lin_colors(i,:))
	else
		histogram(data{i}(:,1),'normalization','pdf','BinMethod','fd','BinWidth',2,'EdgeColor','none','Displayname',cond_list{i},'FaceColor',lin_colors(i,:))
	end
	pd{i} = fitdist(data{i}(:,1),'GeneralizedExtremeValue');
	xl = xlim; yl = ylim;
	%xval = xl(1):diff(xl)/1000:xl(2); 
	xval = 0:0.25:150; 
	yval = pdf(pd{i},xval);
	hold on;
	lin = plot(xval,yval,'LineWidth',2,'Displayname','Gen Ext Val Fit','Color',lin_colors(i,:));
	set( get( get( lin, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
	%set(gca, 'YScale', 'log')
	%xlabel('Generation times (in mins)')
	ylabel(strjoin({cond_list{i},'PDF'},' '));
	xlim([5 90])
end
%legend('location','northeast','NumColumns',2)
sgtitle('Histogram of Generation times across various media')
%name = 'Taher15_genDist_overlap';
xlabel('Generation time (mins)')
name = 'Taher15_genDist_subplot';
savefig(strcat(name,'.fig'));
print(name,'-dpng','-r600');


% Plot standardized distribution	%% All distributions collapse to one! 
figure;
hold on;
for i = 1:7
	var = data{i}(:,1);
	histogram((var-mean(var))/std(var),'normalization','pdf','BinMethod','fd','BinWidth',0.5,'EdgeColor','none','FaceAlpha',0.5,'Displayname',cond_list{i},'FaceColor',lin_colors(i,:))
	pd_{i} = fitdist((var-mean(var))/std(var),'GeneralizedExtremeValue');
	xl = xlim; yl = ylim;
	xval = xl(1):diff(xl)/1000:xl(2); 
	%xval = 0:0.25:150; 
	yval = pdf(pd_{i},xval);
	hold on;
	lin = plot(xval,yval,'LineWidth',2,'Displayname','Gen Ext Val Fit','Color',lin_colors(i,:));
	set( get( get( lin, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
	%set(gca, 'YScale', 'log')
	xlabel('Generation time (in mins)')
	ylabel('PDF')
end
xlim([-6 10])
title({'Histograms of Standardized Doubling times across various media',''})
legend('location','northeast')

name = 'Taher15_genDist_std';
savefig(strcat(name,'.fig'));
print(name,'-dpng','-r600');


%%========================================
% Analysis of standardized histogram

%% GEV Fit pars vs media richness/growth rate
figure; 
for j = 1:7
	var = data{j}(:,1);
	%var = (var-mean(var))/std(var);
	pd{j} = fitdist(var,'GeneralizedExtremeValue');
	pd_k(j) = pd{j}.k;
	pd_sig(j) = pd{j}.sigma;
	pd_mu(j) = pd{j}.mu;
end

% SUbplot version
f = figure; 
subplot(3,1,1)
plot(x, pd_k','-o','LineWidth',1,'DisplayName','Shape Parameter (k)','Color',lin_colors(1,:))
%ylim([0.05 .25])
xl = xlim; yl = ylim;
%ylim([yl(1)*.75 yl(2)+yl(1)*.25])
ylim([yl(1)+yl(2)*.02 yl(2)*.98])
ylabel({'Shape Parameter','(k)',''})
%legend('location','northeast')

subplot(3,1,2)
plot(x, pd_sig','-o','LineWidth',1,'DisplayName','Scale Parameter (\sigma)','Color',lin_colors(2,:))
%ylim([.55 .75])
xl = xlim; yl = ylim;
%ylim([yl(1)*.75 yl(2)+yl(1)*.25])
ylim([yl(1)*.9 yl(2)+yl(1)*.1])
ylabel({'Scale Parameter','(\sigma)',''})
%legend('location','southeast')

subplot(3,1,3)
plot(x, pd_mu','-o','LineWidth',1,'DisplayName','Location Parameter (\mu)','Color',lin_colors(3,:))
%ylim([-0.6 -.35])
xl = xlim; yl = ylim;
ylim([yl(1)*.75 yl(2)+yl(1)*.25])
%ylim([yl(1)+yl(2)*.02 yl(2)*.98])
ylabel({'Location Parameter','(\mu)',''})
%legend('location','northeast')
sgtitle({'Variation of GEV Distribution fit parameters with richer media',''})

name = 'GEV_vs_richerMedia';
%name = 'GEV_vs_richerMedia_std';
savefig(strcat(name,'.fig'));
print(name,'-dpng','-r600');

%% ====

% Skewness and Kurtosis

%% GEV Fit pars vs media richness/growth rate
for i = 1:7
	var = data{i}(:,1);
	var = (var-mean(var))/std(var);
	sk(i) = skewness(var);
	kurt(i) = kurtosis(var);
end

% SUbplot version
f = figure; 
subplot(2,1,1)
plot(x, pd_k','-o','LineWidth',1,'DisplayName','Skewness','Color',lin_colors(1,:))
%ylim([0.05 .25])
xl = xlim; yl = ylim;
%ylim([yl(1)*.75 yl(2)+yl(1)*.25])
ylim([yl(1)+yl(2)*.02 yl(2)*.98])
ylabel({'Skewness',''})
%legend('location','northeast')

subplot(2,1,2)
plot(x, pd_sig','-o','LineWidth',1,'DisplayName','Kurtosis','Color',lin_colors(2,:))
%ylim([.55 .75])
xl = xlim; yl = ylim;
%ylim([yl(1)*.75 yl(2)+yl(1)*.25])
ylim([yl(1)*.9 yl(2)+yl(1)*.1])
ylabel({'Kurtosis',''})
%legend('location','southeast')

sgtitle({'Variation of Moments of Distribution with richer media',''})

%name = 'SkewKurt_vs_richerMedia';
name = 'SkewKurt_vs_richerMedia_std';
savefig(strcat(name,'.fig'));
print(name,'-dpng','-r600');

%%================================

%% =========================
%% Varying Threshold
%% =========================

% COLORS
lin_colors = lines(20);

figure; hold on;
for j = 1:length(threshold)
	doubDat = x_dat(:,:,i,j)/60;
	[h1,h2] = histcounts(doubDat(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	%plot(h2,h1,'linewidth',2)
	%f = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.3,'EdgeColor',lin_colors(j,:),'LineWidth',4)
	f = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'p =',num2str(p(j))},' ')); % ,'LineWidth',4	
	
%	pd = fitdist(doubDat(:),'GeneralizedExtremeValue');
%	xl = xlim; yl = ylim;
%	xval = xl(1):diff(xl)/1000:xl(2); 
%	yval = pdf(pd,xval);
%	l = plot(xval,yval,'LineWidth',1,'Color',[lin_colors(j,:) 0.5],'Displayname','FrechetDist Fit');
end
xlim([0 200])
ylim([0 0.032])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
title({'Histograms of generation times, for varying p, when Sin = 50k/sec',''})
legend('location','northeast')

savefig(strcat('var_p_Hist_gentimes_Sin',num2str(Sin(i)),'.fig'));
print(strcat('var_p_Hist_gentimes_Sin',num2str(Sin(i))),'-dpng','-r600');

%% =========================
% STANDARDIZED

% i => Sin; j => p, 
i = 5; j = 1;
% COLORS
lin_colors = lines(20);

figure; hold on;
for j = 1:length(p)
	doubDat = x_dat(:,:,i,j)/60;
	doubDat = doubDat(:);
	doubDat = (doubDat - mean(doubDat))/std(doubDat);
	[h1,h2] = histcounts(doubDat(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	%plot(h2,h1,'linewidth',2)
	%f = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.3,'EdgeColor',lin_colors(j,:),'LineWidth',4)
	f = fill(h2,h1,lin_colors(j,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'p =',num2str(p(j))},' ')); % ,'LineWidth',4
	
	
%	pd = fitdist(doubDat(:),'GeneralizedExtremeValue')
%	xl = xlim; yl = ylim;
%	xval = xl(1):diff(xl)/1000:xl(2); 
%	yval = pdf(pd,xval);
%	l = plot(xval,yval,'LineWidth',1,'Color',[lin_colors(j,:) 0.5],'Displayname','FrechetDist Fit');
end
xlim([-2 6])
ylim([0 0.65])
xlabel('Standardized generation time')
ylabel('PDF')
title({'Standardized histograms of generation time, for varying p, when Sin = 50k/sec',''})
legend('location','northeast')

savefig(strcat('var_p_STDHist_gentimes_Sin',num2str(Sin(i)),'.fig'));
print(strcat('var_p_STDHist_gentimes_Sin',num2str(Sin(i))),'-dpng','-r600');
