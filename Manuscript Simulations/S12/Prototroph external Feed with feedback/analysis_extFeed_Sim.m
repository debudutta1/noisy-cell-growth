%% Analysis EXT FEED PROTOTROPH

sec_or_min = 0; % sec
n = 3;
p = 3; % [2,3,5];
n_ext_range = [1:p];
seed_range = 0:4;
%feed_rate = [0:2500:10000];
feed_rate = [0 625 1250 2500:2500:15000];

t_ON = 6; 
%t_OFF = 9;
t_OFF = 37;
threshold = 2e7;
same_thresh = 1;

addpath('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\Auxotroph External feed');

%Sin = 1e3*[5 7.5 11 20 50];
Sin = 1e3*[7.5 10 13 20 50];
%Sin = 1e3*[7.5 11 13 20];

gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
createRec = 0;

base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_extFeed\sizer_toff_37\';
%base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_extFeed\sizer_toff_9';

for j = 1:length(feed_rate)
	for n_ext = n_ext_range
		% If no feed, n_ext doesn't matter. Skip redundant simulations
		if feed_rate(j) == 0 & n_ext > 1
			break;
		end
		for i = 1:length(Sin)
			for k = seed_range+1
				load(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k-1),'.mat'),'div_durs','div_durs_exp','size_bir','size_div','size_bir_exp','size_div_exp');

				div_durs_compiled(:,k,j,n_ext,i) = div_durs_exp;
				size_bir_compiled(:,:,k,j,n_ext,i) = size_bir_exp;
				size_div_compiled(:,:,k,j,n_ext,i) = size_div_exp;
				
			end
		end
	end
end

% Obtain Growth rate of cells from dataset
parpool(5)
growth_rate = nan(length(seed_range),length(feed_rate),length(n_ext_range),length(Sin));
reps = 100;
for j = 1:length(feed_rate)
	for n_ext = n_ext_range
		if j == 1 && n_ext > 1
			break;
		end
		for i = 1:length(Sin)
			parfor k = seed_range+1
				tic; growth_rate(k,j,n_ext,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,n_ext,i)); toc;
			end
		end
	end
end
%save(strcat(base_dir,'var_extFeed_growR.mat'),'growth_rate');

%save(strcat(base_dir,'var_extFeed_fig4.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled');
save(strcat(base_dir,'var_extFeed_fig4.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate');

%% ============================================

%% Load Data
%load(strcat(base_dir,'var_extFeed_Dat.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled');
%load(strcat(base_dir,'var_extFeed_growR.mat'),'growth_rate');
%% ============================================

%% ============================================

%% ============================================
% Plot HEATMAP

% Mean Inter-division time based HEATMAP
for i = 1:length(Sin)
	map_dat = permute(mean(mean(div_durs_compiled(:,:,:,:,i))),[4 3 1 2])/60;
	map_dat(map_dat == 0) = nan;
	f = figure; f.Position(4) = 300;
	
	% Relative Scale
%	map_dat = map_dat/map_dat(1,1);		
	
	heatmap(categorical(string(feed_rate)),categorical(string(n_ext_range)), round(map_dat,1),'colormap',viridis(15))
%	f.Children.CellLabelFormat = '%0.3g';	% Relative Scale
	
	xlabel('Rate of metabolite feed (per second)')%,'FontSize',15)
	ylabel('Number of metabolites fed')%,'FontSize',15)
	title({strjoin({'Mean Inter-division times when Sin =', num2str(Sin(i)),'(per sec)'},' '),''});
	print(strcat(base_dir,'heatmap_p_Sin_',num2str(Sin(i))),'-r600','-dpng');
end

% Growth rate HeatMap
% Mean GRowth rate time based HEATMAP
for i = 1:length(Sin)
	map_dat = permute(mean(growth_rate(:,:,:,i)),[3 2 1]);
	map_dat(map_dat == -Inf) = nan;
	f = figure; f.Position(4) = 300;
	
	% Relative Scale
%	map_dat = map_dat/map_dat(1,1);	
	
	heatmap(categorical(string(feed_rate)),categorical(string(n_ext_range)), round(map_dat,2),'colormap',viridis(15))
%	f.Children.CellLabelFormat = '%0.3g';	% Relative Scale
	f.Children.CellLabelFormat = '%0.3g';
	
	xlabel('Rate of metabolite feed (per second)')%,'FontSize',15)
	ylabel('Number of metabolites fed')%,'FontSize',15)
	title({strjoin({'Mean Growth rate when Sin =', num2str(Sin(i)),'(per sec)'},' '),''});
	print(strcat(base_dir,'growthR_heatmap_p_Sin_',num2str(Sin(i))),'-r600','-dpng');
end
%% ============================================
%% Plot HEATMAP with n_ext = 1

% Mean Inter-division time
n_ext = 1;
map_dat = permute(mean(mean(div_durs_compiled(:,:,:,n_ext,:))),[5 3 1 2 4])/60;
f = figure; %f.Position(4) = 340;
heatmap(categorical(string(feed_rate)),categorical(string(Sin)), round(map_dat,1),'colormap',viridis(15))
xlabel('Rate of metabolite feed (per second)')%,'FontSize',15)
ylabel('Rate of substrate import per cascade (per sec)')%,'FontSize',15)
title({'Heatmap of mean Inter-division times when one metabolite is fed externally',''});
print(strcat(base_dir,'heatmap_p_Sin_feed1'),'-r600','-dpng');

% Mean Growth rate
n_ext = 1;
map_dat = permute(mean(growth_rate(:,:,n_ext,:)),[4 2 1 3]);
f = figure; %f.Position(4) = 340;
heatmap(categorical(string(feed_rate)),categorical(string(Sin)), round(map_dat,2),'colormap',viridis(15))
f.Children.CellLabelFormat = '%0.3g'
xlabel('Rate of metabolite feed (per second)')%,'FontSize',15)
ylabel('Rate of substrate import per cascade (per sec)')%,'FontSize',15)
title({'Heatmap of mean Growth rate when one metabolite is fed externally',''});
print(strcat(base_dir,'growthR_heatmap_p_Sin_feed1'),'-r600','-dpng');

%% ============================================
%% Plot HeatMap with histograms

%% Tight SubPlot Specifications
a1 = 0; % a1 - controls the horizontal space between subplots
a2 = 0; % a2 - controls the vertical space between subplots
b1 = 0.1; % b1 - controls the space at the bottom of the entire subplot group
b2 = 0.1; % b2 - controls the space at the top of the entire subplot group
c1 = 0.1; % c1 - controls the space at the left of the entire subplot group
c2 = 0.1; % c2 - controls the space at the right of the entire subplot group
A = [a1 a2];	
B = [b1 b2];
C = [c1 c2];

%%% Feed rate vs Sin
f = figure;
f.Position(4) = 480;
[sp, pos] = tight_subplot(length(feed_rate),length(Sin),A,B,C);
edges = 0:0.025:5;
n_ext = 1;
% Define Colours
cdat = permute(mean(mean(div_durs_compiled(:,:,:,n_ext,:))),[5 3 1 2 4])/60;
cmin = min(cdat(:));
cmax = max(cdat(:));
cmap = viridis(15);
index = fix((cdat-cmin)/(cmax-cmin)*(length(cmap)-1))+1; %A
RGB1 = ind2rgb(index,cmap);
% Permute to get 3D array in proper orientation
RGB2 = permute(RGB1, [2 3 1]);

count = 0;
for j = 1:length(feed_rate)
	for i = 1:length(Sin)
	
		count = count + 1;
			
		[h1,h2] = histcounts(div_durs_compiled(:,:,j,n_ext,i)/3600,edges,'normalization','pdf');
		h2 = h2(1:end-1)+diff(h2)/2;	% Mean of the edges
			
		% Decide if line colour needs to be lightened or darkened, and use function to do the same
		darken = 0;
		if mean(RGB2(j,:,i)) >= 0.45 
			darken = 1;
		end
		axes(sp(count));
		lin_colour = darken_lighten(RGB2(j,:,i),darken, 0.75);
		plot(h2,h1,'linewidth',2,'color',lin_colour)
		%histogram(div_durs_compiled(:,:,j,l,i)/3600,'normalization','pdf')
		xlim([0 3.5])
		ylim([0 1.05*sp(count).YAxis.Limits(2)])
		sp(count).Color = RGB2(j,:,i);
		text(2.2,0.75*sp(count).YAxis.Limits(2),num2str(mean(mean(div_durs_compiled(:,:,j,n_ext,i)))/60,'%0.1f'),'color', lin_colour,'fontweight','bold')
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
ylabel(han,{'Rate of metabolite feed (per sec)',''});
xlabel(han,'Rate of Substrate import per cascade (per sec)');
title(han,'HeatMap with histograms for Prototroph fed metabolites');
%sgtitle('HeatMap with histograms');

% Use subplot xlabel, ylabel to mark secrete_rate and feed_rate
for i = 1:length(Sin)
	sp(length(sp) - length(Sin) + i).XLabel.String = num2str(Sin(i));
end
for i = 1:length(feed_rate)
	sp((i-1)*length(Sin)+1).YLabel.String = num2str(feed_rate(i));
end
% Add the colourmap
han.Colormap = viridis(15)
for i = 1:11
	c_lbl(i) = string(num2str(cmin + (i-1)*(cmax-cmin)/10,'%.1f'));
end
colorbar(han,'Ticks',0:0.1:1,'TickLabels',c_lbl,'location','manual','position',[.91 0.125 0.025 .75]);

f.InvertHardcopy = 'off';
f.Color = 'White';
print(strcat(base_dir,'Hist_heat_p_Sin_1feed'),'-r600','-dpng');


% Feed rate vs n_ext

for i = 1:length(Sin)
	f = figure;
	f.Position = [400 200 520 465];
	[sp, pos] = tight_subplot(length(feed_rate)-1,length(n_ext_range),A,B,C);
	edges = 0:0.025:5;
	
	% Define Colours
	cdat = permute(mean(mean(div_durs_compiled(:,:,2:length(feed_rate),:,i))),[3 4 1 2])/60;
	%cdat = permute(mean(mean(div_durs_compiled(:,:,2:length(feed_rate),:,i))),[4 3 1 2])/60;
	cmin = min(cdat(:));
	cmax = max(cdat(:));
	cmap = viridis(15);
	index = fix((cdat-cmin)/(cmax-cmin)*(length(cmap)-1))+1; %A
	RGB1 = ind2rgb(index,cmap);
	% Permute to get 3D array in proper orientation
	RGB2 = permute(RGB1, [1 3 2]);

	count = 0;
	for j = 2:length(feed_rate)
		for n_ext = 1:length(n_ext_range)
		
			count = count + 1;
			
			[h1,h2] = histcounts(div_durs_compiled(:,:,j,n_ext,i)/3600,edges,'normalization','pdf');
			h2 = h2(1:end-1)+diff(h2)/2;	% Mean of the edges
			
			% Decide if line colour needs to be lightened or darkened, and use function to do the same
			darken = 0;
			if mean(RGB2(j-1,:,n_ext)) >= 0.45 
			%if mean(reshape(RGB1(j-1,n_ext,:),1,3)) >= 0.45 
				darken = 1;
			end
			axes(sp(count));
			lin_colour = darken_lighten(RGB2(j-1,:,n_ext),darken, 0.75);
			%lin_colour = darken_lighten(reshape(RGB1(j-1,n_ext,:),1,3),darken, 0.75);
			plot(h2,h1,'linewidth',2,'color',lin_colour)
			%histogram(div_durs_compiled(:,:,j,l,i)/3600,'normalization','pdf')
			xlim([0 3])
			sp(count).YLim(2) = sp(count).YLim(2)*1.1;
			%sp(count).Color = reshape(RGB1(j-1,n_ext,:),1,3);
			sp(count).Color = RGB2(j-1,:,n_ext);
			text(2.2,0.75*sp(count).YAxis.Limits(2),num2str(mean(mean(div_durs_compiled(:,:,j,n_ext,i)))/60,'%0.1f'),'color', lin_colour,'fontweight','bold')
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
	ylabel(han,{'Rate of metabolite feed (per second)',''});
	xlabel(han,'Number of metabolites fed');
	title(han,strjoin({'HeatMap with histograms, Sin =',num2str(Sin(i)),'(per sec)'},' '));
	
	% Use subplot xlabel, ylabel to mark secrete_rate and feed_rate
	for i1 = 1:length(n_ext_range)
		sp(length(sp) - length(n_ext_range) + i1).XLabel.String = num2str(n_ext_range(i1));
	end
	for i1 = 2:length(feed_rate)
		sp((i1-2)*length(n_ext_range)+1).YLabel.String = num2str(feed_rate(i1));
	end
	% Add the colourmap
	han.Colormap = viridis(15)
	for i1 = 1:11
		c_lbl(i1) = string(num2str(cmin + (i1-1)*(cmax-cmin)/10,'%.1f'));
	end
	colorbar(han,'Ticks',0:0.1:1,'TickLabels',c_lbl,'location','manual','position',[.91 0.125 0.025 .75]);

	f.InvertHardcopy = 'off';
	f.Color = 'White';

	print(strcat(base_dir,'Hist_heat_feedR_nExt_Sin_',num2str(Sin(i))),'-r600','-dpng');	
end





%% ============================================
% Plot line graph for n_ext = 1, mean inter_div vs Sin, for different feed
n_ext = 1;
f = figure; hold on;
for j = 1:length(feed_rate)
	%f.Position(3:4) = [560 460];
	plot(Sin, reshape(mean(mean(div_durs_compiled(:,:,j,n_ext,:)))/60,length(Sin),[]),'-o','linewidth',2,'displayname',strjoin({'Feed rate','=',num2str(feed_rate(j))},' '));
end
title({'Variation of mean Inter-division time with substrate import rate','and varying metabolite feed rate',''})
xlabel({'Rate of Substrate Import', 'per concurrent cascade (# per second)',''})
ylabel('Mean Inter-division time (mins)')
legend('location','northeast')
%xlim([0.6e4 2.25e4]);
f.Children(2).YGrid = 'on';
print(strcat(base_dir,'interDiv_linePlot_Sin_1feed'),'-r600','-dpng');

%% ============================================
% Plot line graph for n_ext = 1, mean inter_div vs feed, for different Sin
n_ext = 1;
f = figure; hold on;
for i = 1:length(Sin)
	%f.Position(3:4) = [560 460];
	plot(feed_rate, reshape(mean(mean(div_durs_compiled(:,:,:,n_ext,i)))/60,length(feed_rate),[]),'-o','linewidth',2,'displayname',strjoin({'S_{in}','=',num2str(Sin(i))},' '));
end
title({'Variation of mean Inter-division time with metabolite feed rate','and varying substrate import rate',''})
xlabel({'Rate of metabolite feed (per sec)',''})
ylabel('Mean Inter-division time (mins)')
legend('location','northeast','numcolumns',2)
%xlim([0.6e4 2.25e4]);
f.Children(2).YGrid = 'on';
print(strcat(base_dir,'interDiv_linePlot_1feed_Sin'),'-r600','-dpng');


%% ============================================

% Plot line graph for n_ext = 1, mean Growth rate vs Sin, for different feed
n_ext = 1;
f = figure; hold on;
for j = 1:length(feed_rate)
	%f.Position(3:4) = [560 460];
	plot(Sin, reshape(nanmean(growth_rate(:,j,n_ext,:)),length(Sin),[]),'-o','linewidth',2,'displayname',strjoin({'Feed rate','=',num2str(feed_rate(j))},' '));
end
title({'Variation of mean Growth rate with Substrate import rate','and varying metabolite feed rate',''})
xlabel({'Rate of Substrate Import', 'per concurrent cascade (# per second)',''})
ylabel('Mean Growth rate (per hour)')
legend('location','southeast')
%xlim([0.6e4 2.25e4]);
%ylim([45 85])
f.Children(2).YGrid = 'on';
print(strcat(base_dir,'growR_linePlot_Sin_p_1feed'),'-r600','-dpng');

%% ============================================

% Plot line graph for n_ext = 1, mean Growth rate vs feed, for different Sin
n_ext = 1;
f = figure; hold on;
for i = 1:length(Sin)
	%f.Position(3:4) = [560 460];
	plot(feed_rate, reshape(nanmean(growth_rate(:,:,n_ext,i)),length(feed_rate),[]),'-o','linewidth',2,'displayname',strjoin({'S_{in}','=',num2str(Sin(i))},' '));
end
title({'Variation of mean Growth rate with metabolite feed rate','and varying Substrate import rate',''})
xlabel({'Rate of metabolite feed (per sec)',''})
ylabel('Mean Growth rate (per hour)')
legend('location','southeast','numcolumns',2)
%xlim([0.6e4 2.25e4]);
%ylim([0.52 0.87])
f.Children(2).YGrid = 'on';
print(strcat(base_dir,'growR_linePlot_1feed_varSin'),'-r600','-dpng');