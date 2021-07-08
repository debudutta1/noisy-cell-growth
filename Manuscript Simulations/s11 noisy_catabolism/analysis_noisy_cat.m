% analysis noisy

% Parameters
Sin = 1e3*[1.5 5.5 6.5 7 7.5 10 20];
p = [1 2 3];
seed_range = 0:2;

% Compile Data from Raw Simulation data
base_dir = 'D:\Debu Simulations\Sep 2020\noisy_catabolism_n2\';

% SIZER
fold_name = '';
for j = 1:length(p)
	for i = 1:length(Sin)
		for k = seed_range
			clear div_durs div_durs_exp sim_vars

			load(strcat(base_dir,fold_name,'\var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs_exp');
			
			x_dat(k+1,:,i,j) = div_durs_exp;
			
%			x_mean(k+1,i,j) = mean(div_durs_exp);
		end
	end
end

% Calculate growth rate from the inter-div_durs_compiled, by fit to an exponential growth curve
% Obtain Growth rate of cells from dataset
gens = 13;
parpool('local',3);
growth_rate = nan(length(seed_range),length(Sin),length(p));
reps = 100;
for j = 1:length(p)
	for i = 1:length(Sin)
		parfor k = seed_range+1
			tic; growth_rate(k,i,j) = exp_grow_rate(reps, gens, x_dat(k,:,i,j)); toc;
		end
	end
end
save('divD_growR_noisy_cat.mat','growth_rate', 'x_dat','Sin','p','seed_range');

%% =========================

%% load data p = 1 2 3 5 10 limited
org = load('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Monod\divD_growR_vary_sin_p_Sizer.mat');
%% =========================


%% =========================
%% Fig - Histogram Fit with GEV
%% =========================


% COLORS
lin_colors = lines(20);


% Plot Histogram ORG (n = 3) vs Noisy Catabolism (n = 2).
%===============================
fg = figure; hold on;
% Sin = 10000 (i=10); p = 1 (j=1)
i = 6; j = 3;
i = 1; j = 3;
doubDat_1 = x_dat(:,:,i,j)/60;
%doubDat1 = doubDat_1(:)/mean(doubDat_1(:));
% Sin = 10000 (i=10); p = 1 (j=1)
i = 11; j = 3;
i = 1; j = 3;
doubDat_2 = org.x_dat(:,:,i,j)/60;
%doubDat2 = doubDat_2(:)/mean(doubDat_2(:));

h1 = histogram(doubDat_1,'normalization','pdf','FaceAlpha',0.5,'EdgeColor', lin_colors(1,:),'DisplayName',strjoin({'Noisy catabolism (n = 2)','Mean =',num2str(mean(doubDat_1(:)),'%.1f'),'CV =',num2str(std(doubDat1),'%.2g')},' ')); %'DisplayStyle', 'stairs','LineWidth',2

h2 = histogram(doubDat_2,'normalization','pdf','FaceAlpha',0.5,'EdgeColor', lin_colors(2,:),'DisplayName',strjoin({'Original model (n = 3)','Mean =',num2str(mean(doubDat_2(:)),'%.1f'),'CV =',num2str(std(doubDat2),'%.2g')},' ')); %'DisplayStyle', 'stairs','LineWidth',2

%xlim([0 200])
xlabel('Cell generation times (mins)')
ylabel('PDF')
legend('location','northeast')

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = strcat('Hist_noisyCat_n2_vs_org_n3_sin10000');
name = strcat('Hist_noisyCat_n2_vs_org_n3_sin1500');
set(fg, 'Color', 'none');
export_fig(name, '-png', '-r600', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));

% Plot Histogram ORG (n = 3) vs Noisy Catabolism (n = 2) -- RESCALED TO MEAN
%===============================
fg = figure; hold on;
% Sin = 10000 (i=10); p = 1 (j=1)
i = 6; j = 3;
i = 1; j = 3;
doubDat_1 = x_dat(:,:,i,j)/60;
doubDat1 = doubDat_1(:)/mean(doubDat_1(:));
% Sin = 10000 (i=10); p = 1 (j=1)
i = 11; j = 3;
i = 1; j = 3;
doubDat_2 = org.x_dat(:,:,i,j)/60;
doubDat2 = doubDat_2(:)/mean(doubDat_2(:));

h1 = histogram(doubDat1,'normalization','pdf','FaceAlpha',0.5,'EdgeColor', lin_colors(1,:),'DisplayName',strjoin({'Noisy catabolism (n = 2)','Mean =',num2str(mean(doubDat_1(:)),'%.1f'),'CV =',num2str(std(doubDat1),'%.2g')},' ')); %'DisplayStyle', 'stairs','LineWidth',2

h2 = histogram(doubDat2,'normalization','pdf','FaceAlpha',0.5,'EdgeColor', lin_colors(2,:),'DisplayName',strjoin({'Original model (n = 3)','Mean =',num2str(mean(doubDat_2(:)),'%.1f'),'CV =',num2str(std(doubDat2),'%.2g')},' ')); %'DisplayStyle', 'stairs','LineWidth',2

%xlim([0 200])
xlabel('Rescaled generation times')
ylabel('PDF')
legend('location','northeast')

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = strcat('Hist_noisyCat_n2_vs_org_n3_sin10000');
name = strcat('Hist_noisyCat_n2_vs_org_n3_sin1500');
set(fg, 'Color', 'none');
export_fig(name, '-png', '-r600', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));


%===============================
% Monods law plot
%===============================
red_scale= 0.75;
%Sin = 1e3*[1.5 5.5 6.5 7 7.5 10 20];
%org.Sin = 1e3*[1.5 3 4 5 5.5 6 6.5 7 7.5 9 10 11 13 20 50];
sel = [1 5 7 8 9 11 14];

fg = figure; hold on;
%for p1 = 1:length(p)
for p1 = 1:3
	plot(org.Sin(sel), reshape(mean(org.growth_rate(:,sel,p1),1),length(org.Sin(sel)),[]),'-o','Displayname',strjoin({'Org Model, p =',num2str(org.p(p1))},' '),'linewidth',2,'color',[lin_colors(p1,:) 0.25])
	
	plot(Sin, reshape(mean(growth_rate(:,:,p1),1),length(Sin),[]),'-<','Displayname',strjoin({'Noisy Cat., p =',num2str(p(p1))},' '),'linewidth',2,'color',[lin_colors(p1,:) 1])
	%SEM = std(growth_rate(:,sel,p1),1);
	%errorbar(Sin(sel), reshape(mean(growth_rate(:,sel,p1),1),length(Sin(sel)),[]), SEM,'-o','Displayname',strjoin({'p =',num2str(p(p1))},' '),'linewidth',2);
end
legend('location','southeast','NumColumns',2);
%title({'Variation of growth rate with rate of import of Substrate','and number of cascades',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux per cascade (per sec)')

xlim([4700 Sin(end)])
%xl = xlim; yl = ylim;
%xlim([xl(1)*.5 xl(2)*1.05])
%ylim([yl(1)*.85 yl(2)*1.02])

fg.Position(3:4) = fg.Position(3:4)*red_scale;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = 'noisy_cat_vs_Org_monod';
set(fg, 'Color', 'none');
export_fig(name, '-png', '-r600', '-transparent', '-painters')
figure2pdf(strcat(name,'.pdf'));

