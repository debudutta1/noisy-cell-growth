% NEW Analysis Code - 20 March 2020

% New Adder vs Sizer

p = [1 2 3 5 10];
Sin = 1e3*[5 7.5 11 20 50];
seed_range = 0:4;
gens =13;
n = 3;

% Load DATA
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_p_sin\Synced_thres_2e7\';
%base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_p_sin\Synced\';

dir_nm = {'Sizer','Adder'};

%% Compiled storage variables
div_durs_compiled = nan(2^gens-1, length(seed_range), length(p), length(Sin), length(dir_nm)); % (:,k,j,i,adder)
%size_bir_compiled = cell(length(p));%nan(2^gens-1, p, length(seed_range), length(feed_rate), length(secrete_rate), length(Sin),2); % (:,:,k,j,i,adder)
%size_div_compiled = cell(length(p));%nan(2^gens-1, p, length(seed_range), length(feed_rate), length(secrete_rate), length(Sin),2); % (:,:,k,j,i,adder)

%size_div_com = nan(2^gens-1, p, length(seed_range), length(Sin), length(dir_nm)); % (:,:,k,i,adder)


%mrna_prod_compiled = nan(p-1,n,2^gens-1, length(seed_range), length(feed_rate), length(secrete_rate), length(Sin),2); % (:,:,:,k,j,i,adder)
%prot_prod_compiled = nan(p-1,n,2^gens-1, length(seed_range), length(feed_rate), length(secrete_rate), length(Sin),2); % (:,:,:,k,j,i,adder)

%mrna_beg_compiled = nan(p-1, n, 2^gens-1, length(seed_range), length(feed_rate), length(secrete_rate), length(Sin),2); % (:,:,:,k,j,i,adder)
%mrna_end_compiled = nan(p-1, n, 2^gens-1, length(seed_range), length(feed_rate), length(secrete_rate), length(Sin),2); % (:,:,:,k,j,i,adder)
%prot_beg_compiled = nan(p-1, n, 2^gens-1, length(seed_range), length(feed_rate), length(secrete_rate), length(Sin),2); % (:,:,:,k,j,i,adder)
%prot_end_compiled = nan(p-1, n, 2^gens-1, length(seed_range), length(feed_rate), length(secrete_rate), length(Sin),2); % (:,:,:,k,j,i,adder)

%dir_nm = {'Sizer','Adder'};

% adder = 1=> Sizer, adder = 2=> Adder dataset
for j = 1:length(p)
	for adder = 1:2;
		for i = 1:length(Sin)
			for k = seed_range+1
				%clear div_durs div_durs_exp sim_vars size_bir_exp size_div_exp
				load(strcat(base_dir,dir_nm{adder},'\var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k-1),'.mat'),'div_durs_exp','size_bir_exp','size_div_exp','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced');
				
				div_durs_compiled(:,k,j,i,adder) = div_durs_exp;
				%size_bir_com(:,:,k,j,i,adder) = size_bir_exp;
				%size_div_com(:,:,k,j,i,adder) = size_div_exp;
				size_bir_com(:,:,k,i,adder) = size_bir_exp;
				size_div_com(:,:,k,i,adder) = size_div_exp;
				%size_bir_compiled(:,:,k,j,i,adder) = size_bir_exp;
				%size_div_compiled(:,:,k,j,i,adder) = size_div_exp;
					
%				mrna_prod_compiled(:,:,:,k,j,i,adder) = mrna_produced;
				mrna_prod_com(:,:,:,k,i,adder) = mrna_produced;
%				prot_prod_compiled(:,:,:,k,j,i,adder) = prot_produced;
				prot_prod_com(:,:,:,k,i,adder) = prot_produced;
					
				mrna_beg_com(:,:,:,k,i,adder) = mrna_beg;
				mrna_end_com(:,:,:,k,i,adder) = mrna_end;
				prot_beg_com(:,:,:,k,i,adder) = prot_beg;
				prot_end_com(:,:,:,k,i,adder) = prot_end;
			end
		end
	end
	size_bir_compiled{j} = size_bir_com;
	size_div_compiled{j} = size_div_com;
	mrna_prod_compiled{j} = mrna_prod_com;
	prot_prod_compiled{j} = prot_prod_com;
	mrna_beg_compiled{j} = mrna_beg_com;
	mrna_end_compiled{j} = mrna_end_com;
	prot_beg_compiled{j} = prot_beg_com;
	prot_end_compiled{j} = prot_end_com;
		
	clear size_bir_com size_div_com mrna_prod_com prot_prod_com mrna_beg_com mrna_end_com prot_beg_com prot_end_com;
end
%save(strcat(base_dir,'sin_p_var_2e7thres_divSiz.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled');
save(strcat(base_dir,'sin_p_var_full_2e7thres.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','mrna_beg_compiled','mrna_end_compiled','mrna_prod_compiled','prot_beg_compiled','prot_end_compiled','prot_prod_compiled');

%========================================
% Load data
p = [1 2 3 5 10];
Sin = 1e3*[5 7.5 11 20 50];
seed_range = 0:4;
gens =13;
n = 3;
dir_nm = {'Sizer','Adder'};

base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_p_sin\Synced_thres_2e7\';
load(strcat(base_dir,'sin_p_var_full_2e7thres.mat'))

%========================================
% Plot mean inter-division times across various Sin and p
for adder = 1:2;
	f = figure; hold on;
	f.Position(3:4) = [560 460];
	for j = 1:length(p)
		plot(Sin, reshape(mean(mean(div_durs_compiled(:,:,j,:,adder)))/60,length(Sin),[]),'-o','linewidth',2,'displayname',strjoin({'p','=',num2str(p(j))},' '));
	end
	title({'Variation of mean Inter-division time with substrate import rate',strjoin({'and number of concurrent cascades for',dir_nm{adder}},' '),''})
	xlabel({'Rate of Substrate Import', 'per concurrent cascade (# per second)',''})
	ylabel('Mean Inter-division time (mins)')
	legend('location','northeast')
	xlim([0.25e4 5.25e4]);
	f.Children(2).YGrid = 'on';
	
	%% print
	print(strcat(base_dir,'interDiv_linePlot_Sin_p_',dir_nm{adder}),'-r300','-dpng');
end

%========================================
% Combined Plot with both Adder and Sizer for mean inter-division times across various Sin and p
f = figure; hold on;
adder = 1;
f.Position(3:4) = [560 460];
for j = 1:length(p)
	plot(Sin, reshape(mean(mean(div_durs_compiled(:,:,j,:,adder)))/60,length(Sin),[]),'-o','linewidth',2,'displayname',strjoin({'p','=',num2str(p(j))},' '));
	lin_colour(j,:) = f.Children.Children(1).Color;
end
adder = 2;
for j = 1:length(p)
	plot(Sin, reshape(mean(mean(div_durs_compiled(:,:,j,:,adder)))/60,length(Sin),[]),'--o','linewidth',2,'displayname',strjoin({'p','=',num2str(p(j))},' '),'color',lin_colour(j,:));
end

title({'Variation of mean Inter-division time with substrate import rate',strjoin({'and number of concurrent cascades for',dir_nm{adder}},' '),''})
xlabel({'Rate of Substrate Import', 'per concurrent cascade (# per second)',''})
ylabel('Mean Inter-division time (mins)')
legend('location','northeast')
xlim([0.25e4 5.25e4]);
f.Children(2).YGrid = 'on';

print(strcat(base_dir,'interDiv_linePlot_Sin_p_adderVSsize'),'-r300','-dpng');

%========================================
% HeatMap

% Mean Inter-division time based HEATMAP
for adder = 1:2
	figure; heatmap(categorical(string(Sin)),categorical(string(p)),permute(mean(mean(div_durs_compiled(:,:,:,:,adder))),[3 4 1 2])/60,'colormap',viridis(15))
	xlabel('Rate of substrate import per cascade (per second)')%,'FontSize',15)
	ylabel('Number of concurrent cascades in the cell')%,'FontSize',15)
	title({strjoin({'Mean Inter-division times when using', dir_nm{adder},'(in mins)'},' '),''});
	print(strcat(base_dir,'heatmap_p_Sin_',dir_nm{adder}),'-r300','-dpng');
end


%========================================
% Extract table of data points
mean_dat = [];
for adder = 1:2;
	for j = 1:length(p)
		mean_dat(j,:,adder) = reshape(mean(mean(div_durs_compiled(:,:,j,:,adder)))/60,length(Sin),[]);
	end
end
display(mean_dat)

%% ===============================
% Extract colourMap for a heatmap using the data
adder = 1;
cdat = mean_dat(:,:,adder);
cmin = min(cdat(:));
cmax = max(cdat(:));
cmap = viridis(15);
index = fix((cdat-cmin)/(cmax-cmin)*(length(cmap)-1))+1; %A

RGB1 = ind2rgb(index,cmap);
% Permute to get 3D array in proper orientation
RGB2 = permute(RGB1, [1 3 2]);
%========================================
% Plot a mean inter-division time subplot hisogram heatmap

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

%% Plot a mean inter-division time subplot hisogram heatmap
f = figure;
[sp, pos] = tight_subplot(length(p),length(Sin),A,B,C);
edges = 0:0.1:5;
count = 0;
%adder = 1;
for j = 1:length(p)
	for i = 1:length(Sin)
		count = count + 1;
		
		[h1,h2] = histcounts(div_durs_compiled(:,:,j,i,adder)/3600,edges,'normalization','pdf');
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
		xlim([0 5])
		%ylim([0 1])
		sp(count).Color = RGB2(j,:,i);
		text(2.5,0.75*sp(count).YAxis.Limits(2),num2str(mean(mean(div_durs_compiled(:,:,j,i,adder)))/60,'%0.1f'),'color', lin_colour,'fontweight','bold')
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
ylabel(han,{'Number of concurrent cascades in the cell',''});
xlabel(han,'Rate of Substrate import per cascade (per second)');
title(han,strjoin({'HeatMap with histograms for simulation w/',dir_nm{adder}},' '));
%sgtitle('HeatMap with histograms');

% Use subplot xlabel, ylabel to mark secrete_rate and feed_rate
for i = 1:length(Sin)
	sp(length(sp) - length(Sin) + i).XLabel.String = num2str(Sin(i));
end
for i = 1:length(p)
	sp((i-1)*length(Sin)+1).YLabel.String = num2str(p(i));
end
% Add the colourmap
han.Colormap = viridis(15)
for i = 1:11
	c_lbl(i) = string(num2str(cmin + (i-1)*(cmax-cmin)/10,'%.1f'));
end
colorbar(han,'Ticks',0:0.1:1,'TickLabels',c_lbl,'location','manual','position',[.91 0.125 0.025 .75]);

f.InvertHardcopy = 'off';
f.Color = 'White';

for i = 1:length(sp)
	sp(i).InvertHardcopy = 'off';
end

print(strcat(base_dir,'Hist_heat_p_Sin_',dir_nm{adder}),'-r300','-dpng');

%========================================
% Plot tight subplot block for Adder vs Sizer overlapp histograms, without color
%% Plot a mean inter-division time subplot hisogram heatmap
f = figure;
[sp, pos] = tight_subplot(length(p),length(Sin),A,B,C);
edges = 0:0.1:5;
count = 0;
for j = 1:length(p)
	for i = 1:length(Sin)
		count = count + 1;
		for adder = [2 1]
			
			[h1,h2] = histcounts(div_durs_compiled(:,:,j,i,adder)/3600,edges,'normalization','pdf');
			h2 = h2(1:end-1)+diff(h2)/2;	% Mean of the edges
			
			lin_colour = [	0.8500    0.3250    0.0980	
							0    0.4470    0.7410];
			
			axes(sp(count)); hold on;
			plot(h2,h1,'linewidth',2,'color',lin_colour(adder,:))
			%histogram(div_durs_compiled(:,:,j,l,i)/3600,'normalization','pdf')
%			if intersect([16 21 22],count)
				xlim([0 5])
%			else
%				xlim([0 3])
%			end
			
			%ylim([0 1])
			%sp(count).Color = RGB2(j,:,i);
			pos1 = [0.85 0.65];
			text(2.7,pos1(adder)*sp(count).YAxis.Limits(2),num2str(mean(mean(div_durs_compiled(:,:,j,i,adder)))/60,'%0.1f'),'color', lin_colour(adder,:),'fontweight','bold')
			sp(count).XAxis.TickLabels = '';
			sp(count).XAxis.TickLength = [0 0];
			sp(count).YAxis.TickLabels = '';
			sp(count).YAxis.TickLength = [0 0];
			%sp(count).XAxis.Visible = 'off';
			%sp(count).YAxis.Visible = 'off';
		end
	end
end
% Added Common title, xlabel, ylabel
han = axes(f,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,{'Number of concurrent cascades in the cell',''});
xlabel(han,'Rate of Substrate import per cascade (per second)');
title(han,'HeatMap with histograms for Adder vs Sizer simulations');
%sgtitle('HeatMap with histograms');

% Use subplot xlabel, ylabel to mark secrete_rate and feed_rate
for i = 1:length(Sin)
	sp(length(sp) - length(Sin) + i).XLabel.String = num2str(Sin(i));
end
for i = 1:length(p)
	sp((i-1)*length(Sin)+1).YLabel.String = num2str(p(i));
end
colormap(lin_colour([2 1],:));
colorbar(han,'Ticks',[0 1],'TickLabels',dir_nm([2 1]),'location','manual','position',[.905 0.45 0.025 .1]);

print(strcat(base_dir,'Hist_var_p_Sin_adder_vs_sizer'),'-r300','-dpng');

%========================================%========================================
%% RESCALED To Mean
% Plot tight subplot block for Adder vs Sizer overlapp histograms, without color
%% Plot a mean inter-division time subplot hisogram heatmap
f = figure;
[sp, pos] = tight_subplot(length(p),length(Sin),A,B,C);
edges = 0:0.1:5;
count = 0;
for j = 1:length(p)
	for i = 1:length(Sin)
		count = count + 1;
		for adder = [2 1]
			
			dat = div_durs_compiled(:,:,j,i,adder)/mean(reshape(div_durs_compiled(:,:,j,i,adder),1,[]));
			[h1,h2] = histcounts(dat(:),edges,'normalization','pdf');
			h2 = h2(1:end-1)+diff(h2)/2;	% Mean of the edges
			
			lin_colour = [	0.8500    0.3250    0.0980	
							0    0.4470    0.7410];
			
			axes(sp(count)); hold on;
			plot(h2,h1,'linewidth',2,'color',lin_colour(adder,:))
			%histogram(div_durs_compiled(:,:,j,l,i)/3600,'normalization','pdf')
%			if intersect([16 21 22],count)
				xlim([0 5])
%			else
%				xlim([0 3])
%			end
			
			%ylim([0 1])
			%sp(count).Color = RGB2(j,:,i);
			pos1 = [0.85 0.65];
			text(2.7,pos1(adder)*sp(count).YAxis.Limits(2),num2str(std(dat(:)),'%0.2f'),'color', lin_colour(adder,:),'fontweight','bold')
			sp(count).XAxis.TickLabels = '';
			sp(count).XAxis.TickLength = [0 0];
			sp(count).YAxis.TickLabels = '';
			sp(count).YAxis.TickLength = [0 0];
			%sp(count).XAxis.Visible = 'off';
			%sp(count).YAxis.Visible = 'off';
		end
	end
end
% Added Common title, xlabel, ylabel
han = axes(f,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,{'Number of concurrent cascades in the cell (p)',''});
xlabel(han,'Rate of Substrate import per cascade (S_{in}) (per second)');
title(han,'Histograms rescaled by the mean for Adder vs Sizer simulations');
%sgtitle('HeatMap with histograms');

% Use subplot xlabel, ylabel to mark secrete_rate and feed_rate
for i = 1:length(Sin)
	sp(length(sp) - length(Sin) + i).XLabel.String = num2str(Sin(i));
end
for i = 1:length(p)
	sp((i-1)*length(Sin)+1).YLabel.String = num2str(p(i));
end
colormap(lin_colour([2 1],:));
colorbar(han,'Ticks',[0 1],'TickLabels',dir_nm([2 1]),'location','manual','position',[.905 0.45 0.025 .1]);

%print(strcat(base_dir,'Hist_meanScaled_var_p_Sin_adder_vs_sizer'),'-r600','-dpng');
set(f, 'Color', 'none');
export_fig('Hist_meanScaled_var_p_Sin_adder_vs_sizer', '-pdf', '-png', '-r600', '-transparent', '-painters')

%========================================
% Compare mean cell size proxies

% Plot mean inter-division times across various Sin and p
for adder = 1:2;
	f = figure; hold on;
	%f.Position(3:4) = [560 460];
	for j = 1:length(p)
		%plot(Sin, reshape(mean(mean(mean(size_bir_compiled{j}(:,:,:,:,adder),2))),length(Sin),[]),'-o','linewidth',2,'displayname',strjoin({'p','=',num2str(p(j))},' '));
		%plot(Sin, reshape(mean(mean(mean(size_div_compiled{j}(:,:,:,:,adder),2))),length(Sin),[]),'-o','linewidth',2,'displayname',strjoin({'p','=',num2str(p(j))},' '));
		plot(Sin, reshape(mean(mean(mean(size_div_compiled{j}(:,:,:,:,adder)-size_bir_compiled{j}(:,:,:,:,adder),2))),length(Sin),[]),'-o','linewidth',2,'displayname',strjoin({'p','=',num2str(p(j))},' '));
	end
	%title({'Variation of mean Inter-division time with substrate import rate',strjoin({'and number of concurrent cascades for',dir_nm{adder}},' '),''})
	xlabel({'Rate of Substrate Import', 'per concurrent cascade (# per second)',''})
	ylabel('Mean cell size')
	legend('location','best')
	xlim([0.25e4 5.25e4]);
	f.Children(2).YGrid = 'on';
	
	%% print
	%print(strcat(base_dir,'interDiv_linePlot_Sin_p_',dir_nm{adder}),'-r300','-dpng');
end

%========================================
% Histogram of size comparisons
figure; hold on;
histogram((mean(size_div_compiled{j}(:,:,:,4,adder)-size_bir_compiled{j}(:,:,:,4,adder),2)))
histogram((mean(size_div_compiled{j}(:,:,:,4,adder+1)-size_bir_compiled{j}(:,:,:,4,adder+1),2)))

%========================================
% Scatter plot of inter-division time vs Added size
adder = 1;
j = 3;
figure; 
X1 = div_durs_compiled(:,:,j,4,adder);
%Y1 = mean(size_div_compiled{j}(:,:,:,4,adder)-size_bir_compiled{j}(:,:,:,4,adder),2);
Y1 = mean(size_div_compiled{j}(:,:,:,4,adder),2)-mean(size_bir_compiled{j}(:,:,:,4,adder),2);
%scatter(X1(:),Y1(:));
%scatter(X1(:)/mean(X1(:)),Y1(:)/mean(Y1(:)));
scatter(X1(:)/mean(X1(:)),Y1(:)/2e7);

% inter-division vs log2
adder = 1;
j = 3;
figure; 
X1 = div_durs_compiled(:,:,j,4,adder);
%Y1 = mean(size_div_compiled{j}(:,:,:,4,adder)-size_bir_compiled{j}(:,:,:,4,adder),2);
Y1 = log2(mean(size_div_compiled{j}(:,:,:,4,adder),2)./mean(size_bir_compiled{j}(:,:,:,4,adder),2));
scatter(X1(:),Y1(:));
%scatter(X1(:)/mean(X1(:)),Y1(:)/mean(Y1(:)));

COLOUR BY ADDED SIZE!!!

%========================================
% Demonstrate Implemented Adder and Sizer work
% Plot of p = 1
j = 1;
for adder = 1:2;
	figure; hold on;
	for i = 1:length(Sin)
		%figure;
		X1 = mean(size_bir_compiled{j}(:,:,:,i,adder),2);
		Y1 = mean(size_div_compiled{j}(:,:,:,i,adder),2)-mean(size_bir_compiled{j}(:,:,:,i,adder),2);
		scatter(X1(:),Y1(:));
		%scatter(X1(:)/mean(X1(:)),Y1(:)/mean(Y1(:)));
	end
	ylabel('Added metabolite')
	xlabel('Metabolite level at birth')
	title({strjoin({'Demonstration of cell division mechanism :',dir_nm{adder}},' '),''})
	ylim(xlim);
	print(strcat(base_dir,'demo_divisionMech_',dir_nm{adder}),'-r300','-dpng');
end

%========================================
% Added size vs birth size
for adder = 1:2;
%j = 3;
	for j = 2:length(p)
	figure; hold on;
		for i = 1:length(Sin)
			%figure;
			X1 = mean(size_bir_compiled{j}(:,:,:,i,adder),2);
			Y1 = mean(size_div_compiled{j}(:,:,:,i,adder),2)-mean(size_bir_compiled{j}(:,:,:,i,adder),2);
			%scatter(X1(:),Y1(:));
			scatter(X1(:)/mean(X1(:)),Y1(:)/mean(Y1(:)));
		end
	end
end
%========================================
% Plot birth size vs division size
for adder = 1:2;
%j = 3;
	for j = 1:length(p)
		figure; hold on;
		for i = 1:length(Sin)
			%figure;
			X1 = mean(size_bir_compiled{j}(:,:,:,i,adder),2);
			Y1 = div_durs_compiled(:,:,j,i,adder)/60;
			%scatter(X1(:),Y1(:));
			scatter(X1(:)/mean(X1(:)),Y1(:)/mean(Y1(:)));
		end
		ylabel('Inter-division time (mins)')
		xlabel('Size at birth')
		%title({strjoin({'Demonstration of cell division mechanism :',dir_nm{adder}},' '),''})
	end
end
%========================================
% Plot birth size vs Inter-division time
for adder = 1:2;
%j = 3;
	for j = 1:length(p)
		figure; hold on;
		for i = 1:length(Sin)
			%figure;
			X1 = mean(size_bir_compiled{j}(:,:,:,i,adder),2);
			Y1 = mean(size_div_compiled{j}(:,:,:,i,adder),2);
			%scatter(X1(:),Y1(:));
			scatter(X1(:)/mean(X1(:)),Y1(:)/mean(Y1(:)));
		end
		ylabel('Division Size')
		xlabel('Size at birth')
		%title({strjoin({'Demonstration of cell division mechanism :',dir_nm{adder}},' '),''})
	end
end

%========================================
% Average mrna protein produced

figure; histogram(prot_prod_compiled(1,:))

%========================================
%% Check for constancy of cell size with linear simulations

Extract Linear cell size data and analyze

%========================================
%========================================
%========================================


%========================================
% Obtain Growth rate of cells from dataset
growth_rate = [];
reps = 100;
for adder = 1:length(dir_nm)
	for j = 1:length(p)
		for i = 1:length(Sin)
			for k = seed_range+1
				growth_rate(k,j,i,adder) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,i,adder));
			end
		end
	end
end
%========================================
% Plot mean growth rate across various Sin and p
for adder = 1:2;
	f = figure; hold on;
	f.Position(3:4) = [560 460];
	for j = 1:length(p)
		plot(Sin, reshape(mean(growth_rate(:,j,:,adder)),length(Sin),[]),'-o','linewidth',2,'displayname',strjoin({'p','=',num2str(p(j))},' '));
	end
	title({'Variation of mean Growth rate with substrate import rate',strjoin({'and number of concurrent cascades for',dir_nm{adder}},' '),''})
	xlabel({'Rate of Substrate Import', 'per concurrent cascade (# per second)',''})
	ylabel('Mean Growth rate (per hour)')
	legend('location','southeast')
	xlim([0.25e4 5.25e4]);
	f.Children(2).YGrid = 'on';
	
	print(strcat(base_dir,'growthR_linePlot_Sin_p_',dir_nm{adder}),'-r300','-dpng');
end
%========================================
% Growth rate HeatMap
% Mean GRowth rate time based HEATMAP
for adder = 1:2
	f = figure; 
	heatmap(categorical(string(Sin)),categorical(string(p)),permute(mean(growth_rate(:,:,:,adder)),[2 3 1]),'colormap',viridis(15));
	f.Children.CellLabelFormat = '%0.2g'
	xlabel('Rate of substrate import per cascade (per second)')%,'FontSize',15)
	ylabel('Number of concurrent cascades in the cell')%,'FontSize',15)
	title({strjoin({'Mean Growth rate when using', dir_nm{adder},'(per hour)'},' '),''});
	print(strcat(base_dir,'growthR_heatmap_p_Sin_',dir_nm{adder}),'-r300','-dpng');
end

%========================================
%========================================
%========================================


% Load data and process
sizer = load('sin_p_var_adderDat.mat');
adder = load('sin_p_var_sizerDat.mat');

% Create standardized datasets
for j = 1:length(p)
	for i = 1:length(Sin)
		for k = seed_range+1
			sizer.std_dat{k,i,j} = (sizer.x_dat{k,i,j} - mean(sizer.x_dat{k,i,j}))./std(sizer.x_dat{k,i,j});
			tem_dat = sizer.x_dat{k,i,j}'/60./log2(mean(sizer.x_dsiz{k,i,j},2)./mean(sizer.x_bsiz{k,i,j},2));
			tem_dat = tem_dat(2:end);	% Neglect the 1st data point has size=0!
			sizer.std_doub{k,i,j} = (tem_dat - mean(tem_dat))./std(tem_dat);
			
			adder.std_dat{k,i,j} = (adder.x_dat{k,i,j} - mean(adder.x_dat{k,i,j}))./std(adder.x_dat{k,i,j});
			tem_dat = adder.x_dat{k,i,j}'/60./log2(mean(adder.x_dsiz{k,i,j},2)./mean(adder.x_bsiz{k,i,j},2));
			tem_dat = tem_dat(2:end);	% Neglect the 1st data point has size=0!
			adder.std_doub{k,i,j} = (tem_dat - mean(tem_dat))./std(tem_dat);
		end
	end
end

%========================================
% Plot inter-division time histogram
figure;histogram(x_dat{3,2,2}/60,'normalization','pdf')

% Plot Doubling time histogram
figure;histogram(x_dat{3,2,2}/60./log2(mean(x_dsiz{3,2,2},2)./mean(x_bsiz{3,2,2},2))','normalization','pdf')

%========================================
% Compare Adder vs Sizer

% Histogram of inter-division times
figure; hold on;
histogram(sizer.x_dat{3,2,2}/60,'normalization','pdf','DisplayName','Sizer')
histogram(adder.x_dat{3,2,2}/60,'normalization','pdf','DisplayName','Adder')
xlabel('Cell inter-division times (in mins)')
ylabel('PDF')
title({'Histogram of inter-division times, from Adder vs Sizer model based simulations',''})
legend('location','northeast')
print('InterDiv_dist_ADvsSZ','-dpng','-r300')


figure; hold on;
k = 3; i = 2; j = 2;
tem_dat = sizer.x_dat{k,i,j}'/60./log2(mean(sizer.x_dsiz{k,i,j},2)./mean(sizer.x_bsiz{k,i,j},2));
histogram(tem_dat,'normalization','pdf','DisplayName','Sizer')
tem_dat = adder.x_dat{k,i,j}'/60./log2(mean(adder.x_dsiz{k,i,j},2)./mean(adder.x_bsiz{k,i,j},2));
histogram(tem_dat,'normalization','pdf','DisplayName','Adder')
xlabel('Cell doubling times (in mins)')
ylabel('PDF')
title({'Histogram of cell Doubling times, from Adder vs Sizer model based simulations',''})
legend('location','northeast')
print('Doubling_dist_ADvsSZ','-dpng','-r300')

%========================================

% Plot histogram with varying sin, for each p
k = 1;
for j = 1:length(p)	
	figure; hold on;
	for i = 1:length(Sin)
		%histogram(x_dat{1,i,j}/60,'normalization','pdf');
		%histogram(x_dat{k,i,j}/60./log2(mean(x_dsiz{k,i,j},2)./mean(x_bsiz{k,i,j},2))','normalization','pdf')
		
		%histogram(sizer.x_dat{k,i,j}/60./log2(mean(sizer.x_dsiz{k,i,j},2)./mean(sizer.x_bsiz{k,i,j},2))','normalization','pdf')
		%histogram(sizer.std_doub{k,i,j},'normalization','pdf')
		
		%histogram(adder.x_dat{k,i,j}/60./log2(mean(adder.x_dsiz{k,i,j},2)./mean(adder.x_bsiz{k,i,j},2))','normalization','pdf')
		histogram(adder.std_doub{k,i,j},'normalization','pdf')
	end
end

% Plot histogram with varying p, for each Sin
for i = 1:length(Sin)
	figure; hold on;
	for j = 1:length(p)	
		%histogram(x_dat{1,i,j}/60,'normalization','pdf');
		histogram(x_dat{k,i,j}/60./log2(mean(x_dsiz{k,i,j},2)./mean(x_bsiz{k,i,j},2))','normalization','pdf')
	end
end