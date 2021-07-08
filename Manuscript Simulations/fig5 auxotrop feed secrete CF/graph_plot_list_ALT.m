% Graphs create script-- Aux Feed Secrete Simulations

%% PLOTTING HEATMAPs

% 1.A. Mean Inter-division time based HEATMAP
for i = 1:length(Sin)
	f = figure; 
	heatmap(categorical(string(secRatio)),categorical(string(feed_rate)),permute(nanmean(nanmean(div_durs_compiled(:,:,:,:,i))),[3 4 1 2])/60,'colormap',viridis(15))
	f.Children.CellLabelFormat = '%0.3g';
	xlabel('Secretion Ratio')%,'FontSize',15)
	ylabel('Metabolite Feed Rate (per second)')%,'FontSize',15)
	title({strjoin({'Mean Inter-division times when Sin','=',num2str(Sin(i))},' '),''});
	print(strcat(base_dir,'heatmap_genT_Sin_',num2str(Sin(i))),'-r600','-dpng')
end

% 1.B. PLOT IN RELATIVE SCALE
% Load ExtFeed Prototroph data for relative scale plot
%extfeed = load('D:\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\ExtFeed Prototroph New\var_extFeed_divdur_growR.mat','div_durs_compiled');
extfeed = load('C:\Users\Dibyendu_Lab\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\ExtFeed Prototroph New\var_extFeed_divdur_growR.mat');
[C1, ia1, ib1] = intersect(feed_rate,extfeed.feed_rate);
[C2, ia2, ib2] = intersect(Sin,extfeed.Sin)
% Mean Inter-division time based HEATMAP
for i = 1:length(Sin)
	f = figure; 
	%heatmap(categorical(string(secRatio)),categorical(string(feed_rate)),round(permute(nanmean(nanmean(div_durs_compiled(:,:,:,:,i)))./nanmean(nanmean(extfeed.div_durs_compiled(:,:,end,1,i))),[3 4 1 2]),2),'colormap',viridis(15))
	
	X_labs = (string(reshape(round(mean(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i))))),1,[]))) + strcat(" ",string(char(177))," ") + (string(reshape(round(std(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i))))),1,[])))
	heatmap(categorical(X_labs),categorical(string(feed_rate)),round(permute(nanmean(nanmean(div_durs_compiled(:,:,:,:,ia2(i)))),[3 4 1 2])./mean(extfeed.mean_div_durs(:,ib1,1,ib2(i)))',2),'colormap',viridis(15))
	f.Children.CellLabelFormat = '%0.3g';
	xlabel('Corresponding mean Secretion Rate (per sec)')
	ylabel('Metabolite Feed Rate (per second)')
	title({strjoin({'Relative mean Inter-division times when Sin','=',num2str(Sin(ia2(i)))},' '),''});
	print(strcat(base_dir,'heatmap_Rel_genT_Sin_',num2str(Sin(ia2(i)))),'-r600','-dpng')
end

%% ===============================
close all
%% ===============================

% Calculate growth rate from the inter-div_durs_compiled, by fit to an exponential growth curve
% Obtain Growth rate of cells from dataset
if ~isfile(strcat(base_dir,'auxFeedSecrete_growR.mat'))
	parpool('local',4);
	growth_rate = nan(length(seed_range),length(feed_rate),length(secRatio),length(Sin));
	reps = 100;
	for j = 1:length(feed_rate)
		for l = 1:length(secRatio)
			for i = 1:length(Sin)
				parfor k = seed_range+1
					tic; growth_rate(k,j,l,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,l,i)); toc;
				end
			end
		end
	end
	save(strcat(base_dir,'auxFeedSecrete_growR.mat'),'growth_rate');
else
	load(strcat(base_dir,'auxFeedSecrete_growR.mat'),'growth_rate');
end
%% ======================
% 2.A. Mean Growth rate based HEATMAP
for i = 1:length(Sin)
	f = figure; 
	heatmap(categorical(string(secRatio)),categorical(string(feed_rate)),round(permute(mean(growth_rate(:,:,:,i)),[2 3 1]),2),'colormap',viridis(15));
	f.Children.CellLabelFormat = '%0.3g';
	xlabel('Secretion Ratio')%,'FontSize',15)
	ylabel('Metabolite Feed Rate (per second)')%,'FontSize',15)
	title({strjoin({'Mean Growth Rate when Sin','=',num2str(Sin(i)),'(per sec per cascade)'},' '),''});
	print(strcat(base_dir,'heatmap_growR_Sin_',num2str(Sin(i))),'-r600','-dpng')
end

% 2.B. PLOT RELATIVE SCALE- Relative to the max growth rate for each Sin, when fed at max levels
% Load ExtFeed Prototroph data for relative scale plot
%extfeed = load('D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_extFeed\sizer_toff_37\var_extFeed_growR.mat','growth_rate');
extfeed = load('C:\Users\Dibyendu_Lab\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\ExtFeed Prototroph New\var_extFeed_divdur_growR.mat');
[C1, ia1, ib1] = intersect(feed_rate,extfeed.feed_rate);
[C2, ia2, ib2] = intersect(Sin,extfeed.Sin);
% Mean Growth rate based HEATMAP
for i = 1:length(Sin)
	f = figure; 
	X_labs = (string(reshape(round(mean(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i))))),1,[]))) + strcat(" ",string(char(177))," ") + (string(reshape(round(std(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i))))),1,[])))
	heatmap(categorical(X_labs),categorical(string(feed_rate)),round(permute(mean(growth_rate(:,:,:,ia2(i))),[2 3 1])./mean(extfeed.growth_rate(:,ib1,1,ib2(i)))',2),'colormap',viridis(15));
	f.Children.CellLabelFormat = '%0.3g';
	xlabel('Corresponding mean Secretion Rate (per sec)')%,'FontSize',15)
	ylabel('Metabolite Feed Rate (per second)')%,'FontSize',15)
	title({strjoin({'Mean relative Growth rate when Sin','=',num2str(Sin(ia2(i))),'(per sec per cascade)'},' '),''});
	print(strcat(base_dir,'heatmap_Rel_growR_Sin_',num2str(Sin(ia2(i)))),'-r600','-dpng')
end

%% ===============================

%% 3. Heatmap with Histograms with inter-division times, for each Sin

% Parameter specification for Sub-Plot without space of histograms of inter-division times - tight_subplot
a1 = 0; % a1 - controls the horizontal space between subplots
a2 = 0; % a2 - controls the vertical space between subplots
b1 = 0.1; % b1 - controls the space at the bottom of the entire subplot group
b2 = 0.1; % b2 - controls the space at the top of the entire subplot group
c1 = 0.1; % c1 - controls the space at the left of the entire subplot group
c2 = 0.1; % c2 - controls the space at the right of the entire subplot group
A = [a1 a2];	
B = [b1 b2];
C = [c1 c2];

cmap = viridis(15);
for i = 1:length(Sin)

	% Extract colourMap for a heatmap using the data
	cdat = permute(nanmean(nanmean(div_durs_compiled(:,:,:,:,i))),[3 4 1 2])/60;
	cmin = nanmin(cdat(:));
	cmax = nanmax(cdat(:));
	index = fix((cdat-cmin)/(cmax-cmin)*(length(cmap)-1))+1; %A

	RGB1 = ind2rgb(index,cmap);
	% Permute to get 3D array in proper orientation
	RGB2 = permute(RGB1, [1 3 2]);

	% Plotting the figure and subplots
	f = figure;
	[sp, pos] = tight_subplot(length(feed_rate),length(secRatio),A,B,C);
	count = 0;
	for j = 1:length(feed_rate)
		for l = 1:length(secRatio)	
			count = count + 1;
			
			edges = 0:0.025:4;
%			if i==1 && l >= 4 
%				edges = 0:1:80;
%			end
			
			[h1,h2] = histcounts(div_durs_compiled(:,:,j,l,i)/3600,edges,'normalization','pdf');
			h2 = h2(1:end-1)+diff(h2)/2;	% Mean of the edges
			
			% Decide if line colour needs to be lightened or darkened, and use function to do the same
			darken = 0;
			if mean(RGB2(j,:,l)) >= 0.45 
				darken = 1;
			end
			axes(sp(count));
			lin_colour = darken_lighten(RGB2(j,:,l),darken, 0.75);
			plot(h2,h1,'linewidth',2,'color',lin_colour)
			%histogram(div_durs_compiled(:,:,j,l,i)/3600,'normalization','pdf')
			if i > 1 && l < 4 
				xlim([0 3])
				%txt_xpos = 1.6;
			end
			ylim([0 1.05*max(h1)])
			sp(count).Color = RGB2(j,:,l);
			text(0.533*sp(count).XAxis.Limits(2),0.75*sp(count).YAxis.Limits(2),num2str(nanmean(nanmean(div_durs_compiled(:,:,j,l,i)))/60,'%0.1f'),'color', lin_colour,'fontweight','bold')
			%text(1.6,0.75*sp(count).YAxis.Limits(2),num2str(mean(div_durs_tbl(j,l,:,i),3)/60,'%0.1f'),'color', lin_colour,'fontweight','bold')
			sp(count).XAxis.TickLabels = '';
			sp(count).XAxis.TickLength = [0 0];
			sp(count).YAxis.TickLabels = '';
			sp(count).YAxis.TickLength = [0 0];
			%sp(count).XAxis.Visible = 'off';
			%sp(count).YAxis.Visible = 'off';
		end
	end
	% Added Common title, xlabel, ylabel
	han = axes(f,'visible','off'); 
	han.Title.Visible='on';
	han.XLabel.Visible='on';
	han.YLabel.Visible='on';
	ylabel(han,{'Metabolite Feed Rate (per second)',''});
	xlabel(han,'Secretion Ratio');
	title(han,{'HeatMap with Histogram of Inter-division times (in mins)',strjoin({'for Sin','=',num2str(Sin(i)),'(per sec per cascade)'},' ')});
	%sgtitle('HeatMap with histograms');

	% Use subplot xlabel, ylabel to mark secRatio and feed_rate
	for i1 = 1:length(secRatio)
		sp(length(sp) - length(secRatio) + i1).XLabel.String = num2str(secRatio(i1));
	end
	for i1 = 1:length(feed_rate)
		sp((i1-1)*length(secRatio)+1).YLabel.String = num2str(feed_rate(i1));
	end
	% Add the colourmap
	han.Colormap = viridis(15)
	for i1 = 1:11
		c_lbl(i1) = string(num2str(cmin + (i1-1)*(cmax-cmin)/10,'%.1f'));
	end
	colorbar(han,'Ticks',0:0.1:1,'TickLabels',c_lbl,'location','manual','position',[.91 0.125 0.025 .75]);
	
	f.InvertHardcopy = 'off';
	f.Color = 'White';
%	print(strcat(base_dir,'Hist_heat_feed_Secrete_Sin',num2str(Sin(i))),'-r600','-dpng');
end

%% ===============================
close all
%% ===============================

% 4. Actual mean Secretion rate obtained due to set Secretion Ratio

% Plot absolute mean Secretion rates obtained
f = figure; hold on;
lin_colour = lines(length(Sin));
% Plot for each of the different cases mean secretion levels
for i = 1:length(Sin)
	for j = 1:length(feed_rate)
		%scatter(secRatio, reshape(mean(sec_end_compiled(:,k,j,:,i)./div_durs_compiled(:,k,j,:,i)),1,[]));
		%h(i) = plot(secRatio, reshape(mean(sec_end_compiled(:,k,j,:,i)./div_durs_compiled(:,k,j,:,i)),1,[]),'-o','linewidth', 2,'color',lin_colour(i,:),'displayname',strjoin({'S_{in}','=',num2str(Sin(i))},' '));
		h(i) = plot(secRatio, reshape(mean(mean(sec_end_compiled(:,:,j,:,i)./div_durs_compiled(:,:,j,:,i))),1,[]),'-o','linewidth', 2,'color',lin_colour(i,:),'displayname',strjoin({'S_{in}','=',num2str(Sin(i))},' '));
	end
end
xlabel('Secretion Ratio set')
ylabel('Resultant mean Secretion Rate')
title({'Secretion rate based on Secretion Ratio set',''})
legend(h,'location','northwest')
f.Children(2).YGrid = 'on';
print(strcat(base_dir,'secRate_vs_SecRatio'),'-r600','-dpng');


% Plot absolute mean internal Metabolite production rates observed
f = figure; hold on;
lin_colour = lines(length(Sin));
% Plot for each of the different cases mean secretion levels
for i = 1:length(Sin)
	for j = 1:length(feed_rate)
		reshape(mean(mean(sec_end_compiled(:,:,j,:,i)./div_durs_compiled(:,:,j,:,i))),1,[])
		%h(i) = plot(secRatio, reshape(mean(reshape(size_div_compiled(:,2,k,j,:,i)-size_bir_compiled(:,2,k,j,:,i),8191,[])./reshape(div_durs_compiled(:,k,j,:,i),8191,[])),1,[]),'-o','linewidth', 2,'color',lin_colour(i,:),'displayname',strjoin({'S_{in}','=',num2str(Sin(i))},' '));
		h(i) = plot(secRatio, reshape(mean(mean(permute((size_div_compiled(:,2,:,j,:,i)-size_bir_compiled(:,2,:,j,:,i)),[1 3 4 5 2])./div_durs_compiled(:,:,j,:,i))),1,[]),'-o','linewidth', 2,'color',lin_colour(i,:),'displayname',strjoin({'S_{in}','=',num2str(Sin(i))},' '));
	end
end
xlabel('Secretion Ratio set')
ylabel('Resultant mean internal Metabolite production Rate')
title({'internal Metabolite production rate based on Secretion Ratio set',''})
legend(h,'location','northeast')
f.Children(2).YGrid = 'on';
print(strcat(base_dir,'metabRate_vs_SecRatio'),'-r600','-dpng');


% Plot Secretion Rate vs Internal Metabolite production rate
f = figure; hold on;
plot(1e4*[0.4 2],1e4*[0.4 2],'--','color','k')
lin_colour = lines(length(Sin));
% Plot for each of the different cases mean secretion levels
for i = 1:length(Sin)
	for j = 1:length(feed_rate)
		%scatter(secRatio, reshape(mean(sec_end_compiled(:,k,j,:,i)./div_durs_compiled(:,k,j,:,i)),1,[]));
		h(i) = plot(reshape(mean(mean(permute((size_div_compiled(:,2,:,j,:,i)-size_bir_compiled(:,2,:,j,:,i)),[1 3 4 5 2])./div_durs_compiled(:,:,j,:,i))),1,[]),reshape(mean(mean(sec_end_compiled(:,:,j,:,i)./div_durs_compiled(:,:,j,:,i))),1,[]),'-o','linewidth', 2,'color',lin_colour(i,:),'displayname',strjoin({'S_{in}','=',num2str(Sin(i))},' '));
	end
end
xlabel('Internal Metabolite production rate (per sec)')
ylabel('Mean Secretion Rate (per sec)')
title({'Comparison of mean Secretion rate vs Internal Metabolite production rate',''})
legend(h,'location','northeast')
f.Children(2).YGrid = 'on';
%xlim([3.5e3 1.5e4])
%ylim([3.5e3 1.5e4])
print(strcat(base_dir,'secRate_vs_MetabRate_wTrend'),'-r600','-dpng');



%% ===========================

%% 5. Protein produced ratio between unperturbed and overexpressed pathway. Ratio calculated: i) over entire population, ii) per cell.
prot_ratio = [];
prot_ratio_percell = [];
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Total number of proteins produced as measure
				% Ratio calculation over entire cellular population
				prot_ratio(k,l,j,i) = sum(sum(prot_prod_compiled(2,:,:,k,j,l,i),2),3)/sum(sum(prot_prod_compiled(1,:,:,k,j,l,i),2),3);
			end
		end
	end
end
%figure; histogram(prot_ratio)

% ============
% Modified to explore more

% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(prot_ratio(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(prot_ratio(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(prot_ratio(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(prot_ratio(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(prot_ratio(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(prot_ratio(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(prot_ratio(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(prot_ratio(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(prot_ratio(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han, {'Mean protein production ratio', 'b/w overproduced and normal',''})
sgtitle({'Variation of protein produced ratio (overexpressed to normal)',''});

print(strcat(base_dir,'subplot_protProd_population_Ratio_var_secRatio'),'-r600','-dpng')


%% ============================

%% 5. Protein produced ratio between unperturbed and overexpressed pathway. Ratio calculated: ii) per cell.

% Plot histogram of ratio of protein production difference between the normal and the overproduced pathway
%count = 0;
prot_ratio = [];
prot_ratio_percell = [];
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Ratio calculation over each cell
				% Protein
				prot_ratio_percell(:,k,l,j,i) = reshape(mean(prot_prod_compiled(2,:,:,k,j,l,i),2)./mean(prot_prod_compiled(1,:,:,k,j,l,i),2),1,[]);
				prot_ratio(k,l,j,i) = mean(prot_ratio_percell(:,k,l,j,i));				
			end
		end
	end
end
%figure; histogram(prot_ratio)

% ============
% Modified to explore more

% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(prot_ratio(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(prot_ratio(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(prot_ratio(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(prot_ratio(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(prot_ratio(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(prot_ratio(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(prot_ratio(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(prot_ratio(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(prot_ratio(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han, {'Mean Protein production ratio', 'b/w overproduced and normal',''})
sgtitle({'Variation of protein produced ratio (overexpressed to normal)',''});

print(strcat(base_dir,'subplot_protProd_percell_Ratio_var_secRatio'),'-r600','-dpng')

%% ==========================
%% 5*. INVERTED RATIO - Protein produced ratio between unperturbed and overexpressed pathway. Ratio calculated: i) over entire population, ii) per cell.
prot_ratio = [];
prot_ratio_percell = [];
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Total number of proteins produced as measure
				% Ratio calculation over entire cellular population
				prot_ratio(k,l,j,i) = sum(sum(prot_prod_compiled(1,:,:,k,j,l,i),2),3)/sum(sum(prot_prod_compiled(2,:,:,k,j,l,i),2),3);
			end
		end
	end
end
%figure; histogram(prot_ratio)

% ============
% Modified to explore more

% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(prot_ratio(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(prot_ratio(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(prot_ratio(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(prot_ratio(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(prot_ratio(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(prot_ratio(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(prot_ratio(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(prot_ratio(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(prot_ratio(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han, {'Mean protein production ratio', 'b/w normal and overexpressed',''})
sgtitle({'Variation of protein produced ratio (normal to overexpressed)',''});

print(strcat(base_dir,'subplot_protProd_population_RatioINV_var_secRatio'),'-r600','-dpng')


%% ============================

%% 5*. INVERTED RATIO - Protein produced ratio between unperturbed and overexpressed pathway. Ratio calculated: ii) per cell.

% Plot histogram of ratio of protein production difference between the normal and the overproduced pathway
%count = 0;
prot_ratio = [];
prot_ratio_percell = [];
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Ratio calculation over each cell
				% Protein
				prot_ratio_percell(:,k,l,j,i) = reshape(mean(prot_prod_compiled(1,:,:,k,j,l,i),2)./mean(prot_prod_compiled(2,:,:,k,j,l,i),2),1,[]);
				prot_ratio(k,l,j,i) = mean(prot_ratio_percell(:,k,l,j,i));				
			end
		end
	end
end
%figure; histogram(prot_ratio)

% ============
% Modified to explore more

% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(prot_ratio(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(prot_ratio(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(prot_ratio(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(prot_ratio(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(prot_ratio(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(prot_ratio(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(prot_ratio(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(prot_ratio(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(prot_ratio(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han, {'Mean Protein production ratio', 'b/w normal and overproduced',''})
sgtitle({'Variation of protein produced ratio (normal to overexpressed)',''});

print(strcat(base_dir,'subplot_protProd_percell_RatioINV_var_secRatio'),'-r600','-dpng')

%% ==========================

%% 6.i)a) Compare total metabolite production ratio between Overexpressed pathway and normal pathway (Per Cell)
tot_metab_OV = [];%nan(2^gens-1, length(seed_range), length(feed_rate), length(secRatio), length(Sin)); % (:,k,j,l,i)
metab_ratio = [];%nan(2^gens-1, length(seed_range), length(feed_rate), length(secRatio), length(Sin)); % (:,k,j,l,i)
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Calculation for each cell Total internal amount
				tot_metab_OV(:,k,l,j,i) = sec_end_compiled(:,k,j,l,i) + size_div_compiled(:,2,k,j,l,i) - size_bir_compiled(:,2,k,j,l,i);
				metab_ratio_percell(:,k,l,j,i) = tot_metab_OV(:,k,l,j,i)./(size_div_compiled(:,1,k,j,l,i) - size_bir_compiled(:,1,k,j,l,i));
				metab_ratio(k,l,j,i) = mean(metab_ratio_percell(:,k,l,j,i));
			end
		end
	end
end
%figure; histogram(metab_ratio);
mean(metab_ratio(:))

%% ===============================

% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(metab_ratio(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(metab_ratio(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(metab_ratio(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(metab_ratio(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(metab_ratio(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(metab_ratio(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(metab_ratio(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(metab_ratio(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(metab_ratio(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';
% Total metabolite produced
ylabel(han, {'Mean total metabolite produced ratio', 'b/w overproduced and normal',''})
sgtitle({'Variation of total metabolite produced', 'ratio (overexpressed to normal)',''});
print(strcat(base_dir,'subplot_fluxRatioPERCELL_var_secRatio'),'-r600','-dpng')


% ========================
%% 6.i)b) Compare total metabolite production ratio between Overexpressed pathway and normal pathway (POPULATION)
tot_metab_OV = [];%nan(2^gens-1, length(seed_range), length(feed_rate), length(secRatio), length(Sin)); % (:,k,j,l,i)
metab_ratio = [];%nan(2^gens-1, length(seed_range), length(feed_rate), length(secRatio), length(Sin)); % (:,k,j,l,i)
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Calculation for each cell Total internal amount
				tot_metab_OV(:,k,l,j,i) = sec_end_compiled(:,k,j,l,i) + size_div_compiled(:,2,k,j,l,i) - size_bir_compiled(:,2,k,j,l,i);
				metab_ratio(k,l,j,i) = sum(tot_metab_OV(:,k,l,j,i))./sum(size_div_compiled(:,1,k,j,l,i) - size_bir_compiled(:,1,k,j,l,i));
			end
		end
	end
end
%figure; histogram(metab_ratio);
mean(metab_ratio(:))

%% ===============================

% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(metab_ratio(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(metab_ratio(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(metab_ratio(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(metab_ratio(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(metab_ratio(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(metab_ratio(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(metab_ratio(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(metab_ratio(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(metab_ratio(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';
% Total metabolite produced
ylabel(han, {'Mean total metabolite produced ratio', 'b/w overproduced and normal',''})
sgtitle({'Variation of total metabolite produced', 'ratio (overexpressed to normal)',''});
print(strcat(base_dir,'subplot_fluxRatioPOPULATION_var_secRatio'),'-r600','-dpng')

%% ============================

%% 6.ii)a) Compare Internal metabolite production ratio between Overexpressed pathway and normal pathway (PER CELL)
tot_metab_OV = [];
metab_ratio = [];
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Calculation for each cell - only internal amount of overexpressed metabolite
				tot_metab_OV(:,k,l,j,i) = size_div_compiled(:,2,k,j,l,i) - size_bir_compiled(:,2,k,j,l,i);
				metab_ratio_percell(:,k,l,j,i) = tot_metab_OV(:,k,l,j,i)./(size_div_compiled(:,1,k,j,l,i) - size_bir_compiled(:,1,k,j,l,i));
				metab_ratio(k,l,j,i) = mean(metab_ratio_percell(:,k,l,j,i));
			end
		end
	end
end
%figure; histogram(metab_ratio);
mean(metab_ratio(:))

%% ===============================

% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(metab_ratio(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(metab_ratio(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(metab_ratio(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(metab_ratio(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(metab_ratio(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(metab_ratio(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(metab_ratio(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(metab_ratio(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(metab_ratio(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';

% Internal metabolite produced
ylabel(han, {'Mean Internal metabolite produced ratio', 'b/w overproduced and normal',''})
sgtitle({'Variation of Internal metabolite produced','ratio (overexpressed to normal)',''});
print(strcat(base_dir,'subplot_InternalfluxRatio_PERCELL_var_secRatio'),'-r600','-dpng')

%% 6.ii)b) Compare Internal metabolite production ratio between Overexpressed pathway and normal pathway (POPULATION)
tot_metab_OV = [];
metab_ratio = [];
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Calculation for each cell - only internal amount of overexpressed metabolite
				tot_metab_OV(:,k,l,j,i) = size_div_compiled(:,2,k,j,l,i) - size_bir_compiled(:,2,k,j,l,i);
				metab_ratio(k,l,j,i) = sum(tot_metab_OV(:,k,l,j,i))./sum(size_div_compiled(:,1,k,j,l,i) - size_bir_compiled(:,1,k,j,l,i));
			end
		end
	end
end
%figure; histogram(metab_ratio);
mean(metab_ratio(:))

%% ===============================

% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(metab_ratio(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(metab_ratio(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(metab_ratio(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(metab_ratio(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(metab_ratio(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(metab_ratio(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(metab_ratio(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(metab_ratio(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(metab_ratio(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';

% Internal metabolite produced
ylabel(han, {'Mean Internal metabolite produced ratio', 'b/w overproduced and normal',''})
sgtitle({'Variation of Internal metabolite produced','ratio (overexpressed to normal)',''});
print(strcat(base_dir,'subplot_InternalfluxRatio_POPULATION_var_secRatio'),'-r600','-dpng')

%% ====================

%% 7. Plot SUBPLOT to study variation of Inter-division time variation across various parameters
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				mean_div_dur(k,l,j,i) = mean(div_durs_compiled(:,k,j,l,i))/60;
			end
		end
	end
end
% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(mean_div_dur(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(mean_div_dur(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(mean_div_dur(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(mean_div_dur(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(mean_div_dur(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(mean_div_dur(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(mean_div_dur(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(mean_div_dur(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(mean_div_dur(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';
ylabel(han, {'Mean inter-division time',''})
sgtitle({'Variation of inter-division time across parameters',''});
print(strcat(base_dir,'subplot_divdurs_var_secRatio'),'-r600','-dpng')

%% ====================

%% 8. Plot SUBPLOT to study variation of Growth Rate variation across various parameters
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				mean_div_dur(k,l,j,i) = growthR_tbl(j,l,k,i); 
			end
		end
	end
end
% Plot means across each parameter
f = figure;
% Across secRatio
sp(1) = subplot(2,2,1);
hold on;
plot_dat = mean(mean(mean(mean_div_dur(:,:,:,:),1),3),4);
%plot(secRatio, mean(mean(mean(mean_div_dur(:,:,:,:),1),3),4),'-o','linewidth',2,'color','k','displayname','Mean')
plot(secRatio, plot_dat,'linewidth',2,'color','k','displayname','Mean')
for l = 1:length(secRatio)
	scatter(secRatio(l), plot_dat(l),100,'filled','s')
end
%ylim([1.77 2.01])
xlabel('Secretion Ratio')
%for k = seed_range+1
%	plot(secRatio, reshape(mean(mean(mean_div_dur(k,:,:,:),3),4),1,[]),'-o','displayname',strjoin({'Seed','=',num2str(k-1)},' '))
%end
%legend('location','best')
% Across feed_rate
sp(2) = subplot(2,2,2);
hold on;
plot(feed_rate, reshape(mean(mean(mean(mean_div_dur(:,:,:,:),1),2),4),1,[]),'-o','linewidth',2,'color','k','displayname','Mean')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Metabolite Feed Rate')
for l = 1:length(secRatio)
	plot(feed_rate, reshape(mean(mean(mean_div_dur(:,l,:,:),1),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across Sin
sp(3) = subplot(2,2,3);
hold on;
plot(Sin, reshape(mean(mean(mean(mean_div_dur(:,:,:,:),1),2),3),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Substrate Import Rate')
for l = 1:length(secRatio)
	plot(Sin, reshape(mean(mean(mean_div_dur(:,l,:,:),1),3),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Across different Seeds
sp(4) = subplot(2,2,4);
hold on;
plot(seed_range, reshape(mean(mean(mean(mean_div_dur(:,:,:,:),2),3),4),1,[]),'-o','linewidth',2,'color','k')
%ylim([1.8 2])
%ylim([1.77 2.01])
xlabel('Seed values')
for l = 1:length(secRatio)
	plot(seed_range, reshape(mean(mean(mean_div_dur(:,l,:,:),3),4),1,[]),'-o','linewidth',2,'displayname',strjoin({'Sec Rt','=',num2str(secRatio(l))},' '))
end
%legend('location','best')
% Labelling
han = axes(f,'visible','off'); 
%han.Title.Visible='on';
han.YLabel.Visible='on';

ylabel(han, {'Mean Growth rate',''})
sgtitle({'Variation of Growth rate across parameters',''});
print(strcat(base_dir,'subplot_growthR_var_secRatio'),'-r600','-dpng')


%%% =======================
% 2.B.(II) PLOT RELATIVE SCALE- Relative to the max growth rate for each Sin, when fed at max levels
%% Find and mark simulated conditions where there is net surplus with a faster growth rate

% Load ExtFeed Prototroph data for relative scale plot
%extfeed = load('D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_extFeed\sizer_toff_37\var_extFeed_growR.mat','growth_rate');
extfeed = load('C:\Users\Dibyendu_Lab\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Figures\ExtFeed Prototroph New\var_extFeed_divdur_growR.mat');
[C1, ia1, ib1] = intersect(feed_rate,extfeed.feed_rate);
[C2, ia2, ib2] = intersect(Sin,extfeed.Sin);
% Mean Growth rate based HEATMAP
for i = 1:length(Sin)
	f = figure; 
	heat_dat = round(permute(mean(growth_rate(:,:,:,ia2(i))),[2 3 1])./mean(extfeed.growth_rate(:,ib1,1,ib2(i)))',2);
	% Replace Data that has no surplus, and higher growth rate than extFeedProt, with NaN
	select1 = ((reshape(round(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i)))),length(feed_rate),[]) - feed_rate') >= 0);
	select2 = (heat_dat >= 1);
	select = select1.*select2;
	heat_dat(~select) = nan;
	
	X_labs = (string(reshape(round(mean(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i))))),1,[]))) + strcat(" ",string(char(177))," ") + (string(reshape(round(std(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i))))),1,[])));
	
	heatmap(categorical(X_labs),categorical(string(feed_rate)),heat_dat,'colormap',viridis(15));
	f.Children.CellLabelFormat = '%0.3g';
	xlabel('Corresponding mean Secretion Rate (per sec)')%,'FontSize',15)
	ylabel('Metabolite Feed Rate (per second)')%,'FontSize',15)
	title({strjoin({'Mean relative Growth rate when Sin','=',num2str(Sin(ia2(i))),'(per sec per cascade)'},' '),''});
	print(strcat(base_dir,'heatmap_feasible_Rel_growR_Sin_',num2str(Sin(ia2(i)))),'-r600','-dpng')
end

%% PLOT ONLY SURPLUS HEATMAP
% Mean Growth rate based HEATMAP
for i = 1:length(Sin)
	f = figure; 
	heat_dat = round(permute(mean(growth_rate(:,:,:,ia2(i))),[2 3 1])./mean(extfeed.growth_rate(:,ib1,1,ib2(i)))',2);
	% Replace Data that has no surplus, and higher growth rate than extFeedProt, with NaN
	select = ((reshape(round(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i)))),length(feed_rate),[]) - feed_rate') >= 0);
	heat_dat(~select) = nan;
	
	X_labs = (string(reshape(round(mean(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i))))),1,[]))) + strcat(" ",string(char(177))," ") + (string(reshape(round(std(mean(mean(sec_end_compiled(:,:,:,:,i)./div_durs_compiled(:,:,:,:,i))))),1,[])));
	
	heatmap(categorical(X_labs),categorical(string(feed_rate)),heat_dat,'colormap',viridis(15));
	f.Children.CellLabelFormat = '%0.3g';
	xlabel('Corresponding mean Secretion Rate (per sec)')%,'FontSize',15)
	ylabel('Metabolite Feed Rate (per second)')%,'FontSize',15)
	title({strjoin({'Mean relative Growth rate when Sin','=',num2str(Sin(ia2(i))),'(per sec per cascade)'},' '),''});
	print(strcat(base_dir,'heatmap_surplus_Rel_growR_Sin_',num2str(Sin(ia2(i)))),'-r600','-dpng')
end


%f = figure;
%C1 = [0 2 nan 6; 8 nan 12 14; 16 18 20 22];
%im1 = heatmap(C1,'colormap',viridis(10));
%colormap viridis(10);
%im1 = imagesc(C1);
%hold on;
%C2 = [0 1 0 0; 1 1 0 1; 0 0 0 1];
%colormap [0 0 0; 1 1 1];
%im2 = imagesc(C2);
%im2.AlphaData = .5;
%im2.Parent.XAxis.Visible = 0
%im2.Parent.YAxis.Visible = 0


prot_ratio = [];
for j = 1:length(feed_rate)
	for l = 1:length(secRatio)
		for i = 1:length(Sin)
			for k = seed_range+1
				% Total number of proteins produced as measure
				% Ratio calculation over entire cellular population
				prot_ratio(k,l,j,i) = sum(sum(prot_prod_compiled(2,:,:,k,j,l,i),2),3)/sum(sum(prot_prod_compiled(1,:,:,k,j,l,i),2),3);
			end
		end
	end
end

%% PLOT HEATMAP PROT RATIO POPULATION
% Mean Growth rate based HEATMAP
for i = 1:length(Sin)
	f = figure;
	heatmap(categorical(string(secRatio)),categorical(string(feed_rate)),permute(mean(prot_ratio(:,:,:,i)),[3 2 1]),'colormap',viridis(15));
	f.Children.CellLabelFormat = '%0.3g';
	xlabel('Corresponding mean Secretion Ratio')%,'FontSize',15)
	ylabel('Metabolite Feed Rate (per second)')%,'FontSize',15)
	title({strjoin({'Mean Protein produced ratio when Sin','=',num2str(Sin(i)),'(per sec per cascade)'},' '),''});
	print(strcat(base_dir,'heatmap_pop_protRatio_Sin_',num2str(Sin(i))),'-r600','-dpng')
end
close all