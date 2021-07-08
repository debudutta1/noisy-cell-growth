% Manuscript Figure 1

addpath('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\')

% Linear Simulation of 100 cell divs
Sin = 9e3;
p = 2;
%p = 3;
seed_range = 0;

%gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
gens = 1; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;

%createRec = 0;
createRec = 1;	% Detailed full output
%createRec = 2;	% Detailed data not output, but computed to extract important parameters

same_thresh = 1;
if same_thresh == 1
	threshold = 1e+07;
else
	for i = 1:p
		threshold(i) = 1e+07;
	end
end

% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent

for k = seed_range
	for j = 1:length(p)
		for i = length(Sin):-1:1
			
			rng(k);
			y = parallel_growth_sim_v2_1(gens, Sin(i), n,p(j), T, config, sec_or_min, createRec, threshold, same_thresh);
		end
	end
end
%% SImulation ERROR due to missing mrna_Produced for createRec = 1;

proto_cell =  y.proto_cell;

save('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Analysis\Fig1\proto_cell.mat','proto_cell')
%++++=====================

load('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Analysis\Fig1\proto_cell.mat','proto_cell');

% Plot a subplot(4,1) , with burst, mRNA, Protein, and metabolite
figure;
%range = 1:10;
%range = 10:50;
%range = 37:41;
range = 36:38;

% Burst Prof
h1 = subplot(4,1,1);
hold on;

LW = 0.75;	% LineWidth
% Combine data from cell range
start = []; stop = [];
for i = range
	start = [start proto_cell{i}.burst_dat{1,1}.start/60];
	stop = [stop proto_cell{i}.burst_dat{1,1}.stop/60];
end
beg = proto_cell{range(1)}.metab_prof.x(1)/60;
burst = [beg start stop];
c = [0 ones(1,length(start)) zeros(1,length(stop))];
[a,b] = sort(burst);
burst = a;
c = c(b);
p1 = stairs(burst, c,'linewidth',LW);

% Combine data from cell range
start = []; stop = [];
for i = range
	start = [start proto_cell{i}.burst_dat{2,1}.start/60];
	stop = [stop proto_cell{i}.burst_dat{2,1}.stop/60];
end
beg = proto_cell{range(1)}.metab_prof.x(1)/60;
burst = [beg start stop];
c = [0 ones(1,length(start)) zeros(1,length(stop))];
[a,b] = sort(burst);
burst = a;
c = c(b);
p2 = stairs(burst, c,'linewidth',LW);

%h1.XAxis.Visible = 0
%%=========
% mRNA Prof
h2 = subplot(4,1,2);
hold on;

LW2 = 1;	% LineWidth
% Combine data from cell range
data1.t =[];data2.t =[];data1.v =[];data2.v =[];
for i = range
	data1.t = [data1.t proto_cell{i}.mrna_prof{1,1}.t/60];
	data1.v = [data1.v proto_cell{i}.mrna_prof{1,1}.v];
	data2.t = [data2.t proto_cell{i}.mrna_prof{2,1}.t/60];
	data2.v = [data2.v proto_cell{i}.mrna_prof{2,1}.v];
end
p1 = stairs(data1.t, data1.v,'linewidth',LW2);
hold on;
p2 = stairs(data2.t, data2.v,'linewidth',LW2);
%h2.XAxis.Visible = 0;

%%=========
% Prot Prof
h3 = subplot(4,1,3);
hold on;

% Combine data from cell range
data1.t =[];data2.t =[];data1.v =[];data2.v =[];
for i = range
	data1.t = [data1.t proto_cell{i}.prot_prof{1,1}.t/60];
	data1.v = [data1.v proto_cell{i}.prot_prof{1,1}.v];
	data2.t = [data2.t proto_cell{i}.prot_prof{2,1}.t/60];
	data2.v = [data2.v proto_cell{i}.prot_prof{2,1}.v];
end
p1 = stairs(data1.t, data1.v,'linewidth',LW2);
hold on;
p2 = stairs(data2.t, data2.v,'linewidth',LW2);
%h3.XAxis.Visible = 0;

%%=========
% Metab Prof
h4 = subplot(4,1,4);
hold on;

% Combine data from cell range
data.x =[]; data.y =[];
for i = range
	data.x = [data.x proto_cell{i}.metab_prof.x/60];
	data.y = [data.y proto_cell{i}.metab_prof.y(end-p+1:end,:)];
end
p1 = plot(data.x, data.y','linewidth',LW2);
%h4.XAxis.Visible = 0;

%% Formatting the plot

% Hide Time axis keep one
h1.XAxis.Visible = 0;
h2.XAxis.Visible = 0;
h3.XAxis.Visible = 0;

%x_limits = [2060 2293];
%x_limits = h4.YAxis.Limits;
x_limits = [1166 1260];
h1.XLim = x_limits;
h2.XLim = x_limits;
h3.XLim = x_limits;
h4.XLim = x_limits;

%h4.YLim = [1.5e7 5.5e7];

h4.YLabel = '';

%++++=====================


%++++=====================
%% Alternative Individual Graphs
%++++=====================
% Plot burst, mRNA, Protein, and metabolite
range = 36:38;
lin_weight = 1.5;
%x_limits = [2060 2293];
x_limits = [1166 1260];
threshold = 1e7;

% Burst Prof
h1 = figure;
hold on;

% Combine data from cell range
start = []; stop = [];
for i = range
	start = [start proto_cell{i}.burst_dat{1,1}.start/60];
	stop = [stop proto_cell{i}.burst_dat{1,1}.stop/60];
end
beg = proto_cell{range(1)}.metab_prof.x(1)/60;
burst = [beg start stop];
c = [0 ones(1,length(start)) zeros(1,length(stop))];
[a,b] = sort(burst);
burst = a;
c = c(b);
p1 = stairs(burst, c,'linewidth',lin_weight);

% Combine data from cell range
start = []; stop = [];
for i = range
	start = [start proto_cell{i}.burst_dat{2,1}.start/60];
	stop = [stop proto_cell{i}.burst_dat{2,1}.stop/60];
end
beg = proto_cell{range(1)}.metab_prof.x(1)/60;
burst = [beg start stop];
c = [0 ones(1,length(start)) zeros(1,length(stop))];
[a,b] = sort(burst);
burst = a;
c = c(b);
p2 = stairs(burst, c,'linewidth',lin_weight);

h1.Children.XAxis.Visible = 0;
%h1.Position(3:4) = [1150 80];
h1.Position = [96   400   740    80];
h1.Children.YAxis.FontSize = 11;
h1.Children.YAxis.TickValues = [0 1];
h1.Children.XLim = x_limits;

% export_fig
set(h1, 'Color', 'none');
export_fig burst -pdf -png -r600 -transparent


%%=========
% mRNA Prof
h2 =  figure;
hold on;

% Combine data from cell range
data1.t =[];data2.t =[];data1.v =[];data2.v =[];
for i = range
	data1.t = [data1.t proto_cell{i}.mrna_prof{1,1}.t/60];
	data1.v = [data1.v proto_cell{i}.mrna_prof{1,1}.v];
	data2.t = [data2.t proto_cell{i}.mrna_prof{2,1}.t/60];
	data2.v = [data2.v proto_cell{i}.mrna_prof{2,1}.v];
end
p1 = stairs(data1.t, data1.v,'linewidth',lin_weight);
hold on;
p2 = stairs(data2.t, data2.v,'linewidth',lin_weight);

h2.Position = [96   400   740    100];
h2.Children.YAxis.FontSize = 11;
h2.Children.XLim = x_limits;
h2.Children.XAxis.Visible = 0;
% export_fig
set(h2, 'Color', 'none');
export_fig mRNA -pdf -png -r600 -transparent

%%=========
% Prot Prof
h3 = figure;
hold on;

% Combine data from cell range
data1.t =[];data2.t =[];data1.v =[];data2.v =[];
for i = range
	data1.t = [data1.t proto_cell{i}.prot_prof{1,1}.t/60];
	data1.v = [data1.v proto_cell{i}.prot_prof{1,1}.v];
	data2.t = [data2.t proto_cell{i}.prot_prof{2,1}.t/60];
	data2.v = [data2.v proto_cell{i}.prot_prof{2,1}.v];
end
p1 = stairs(data1.t, data1.v,'linewidth',lin_weight);
hold on;
p2 = stairs(data2.t, data2.v,'linewidth',lin_weight);

h3.Position = [96   400   740    100];
h3.Children.YAxis.FontSize = 11;
h3.Children.XLim = x_limits;
h3.Children.XAxis.Visible = 0;
% export_fig
set(h3, 'Color', 'none');
export_fig Protein -pdf -png -r600 -transparent

%%=========
% Metab Prof
h4 = figure;
hold on;

% Combine data from cell range
data.x =[]; data.y =[];
for i = range
	data.x = [data.x proto_cell{i}.metab_prof.x/60];
	data.y = [data.y proto_cell{i}.metab_prof.y(end-p+1:end,:)];
end
p1 = plot(data.x, data.y','linewidth',lin_weight);
p2 = plot(x_limits, [threshold threshold]*2,'--','linewidth',0.5,'color','k');

h4.Position = [96   400   740    125];
h4.Children.YAxis.FontSize = 11;
h4.Children.XAxis.FontSize = 11;
h4.Children.XLim = x_limits;
h4.Children.YLim = 1e7*[0.8 2.5];

h4.Children.XAxis.TickLabels = cellstr(string(h4.Children.XAxis.TickValues - x_limits(1)))';

% Reorder Children
h4.Children.Children = [ h4.Children.Children(2) h4.Children.Children(3) h4.Children.Children(1) ];

% export_fig
set(h4, 'Color', 'none');
export_fig Metabolite -pdf -png -r600 -transparent


%% ========= BURST PROFILE AS SEPARATE GRAPHS
% Line COLORS
lin_colors = lines(2);

% Burst Prof
h1 = figure;
hold on;

% Combine data from cell range
start = []; stop = [];
for i = range
	start = [start proto_cell{i}.burst_dat{1,1}.start/60];
	stop = [stop proto_cell{i}.burst_dat{1,1}.stop/60];
end
beg = proto_cell{range(1)}.metab_prof.x(1)/60;
burst = [beg start stop];
c = [0 ones(1,length(start)) zeros(1,length(stop))];
[a,b] = sort(burst);
burst = a;
c = c(b);
subplot(2,1,1)
p2 = stairs(burst, c,'linewidth',lin_weight,'color',lin_colors(1,:));

% Combine data from cell range
start = []; stop = [];
for i = range
	start = [start proto_cell{i}.burst_dat{2,1}.start/60];
	stop = [stop proto_cell{i}.burst_dat{2,1}.stop/60];
end
beg = proto_cell{range(1)}.metab_prof.x(1)/60;
burst = [beg start stop];
c = [0 ones(1,length(start)) zeros(1,length(stop))];
[a,b] = sort(burst);
burst = a;
c = c(b);
subplot(2,1,2)
p2 = stairs(burst, c,'linewidth',lin_weight,'color',lin_colors(2,:));

h1.Children(1).XAxis.Visible = 0;
h1.Children(2).XAxis.Visible = 0;
h1.Position = [96   400   740    100];
h1.Children(1).YAxis.FontSize = 11;
h1.Children(2).YAxis.FontSize = 11;
h1.Children(1).YAxis.TickValues = [0 1];
h1.Children(2).YAxis.TickValues = [0 1];
h1.Children(1).XLim = x_limits;
h1.Children(2).XLim = x_limits;

% export_fig
set(h1, 'Color', 'none');
export_fig burst -pdf -png -r600 -transparent




















%++++=====================
% ROUGH
p = 2;
%++++=====================
% Transcription Burst plot
%++++=====================

burst = [proto_cell{1}.burst_dat{1,1}.start/60 proto_cell{1}.burst_dat{1,1}.stop/60];
c = [ones(1,length(proto_cell{1}.burst_dat{1,1}.start)) zeros(1,length(proto_cell{1}.burst_dat{1,1}.stop))];
[a,b] = sort(burst);
burst = a;
c = c(b);

f = figure;
p1 = stairs(burst, c,'linewidth',2);

burst = [proto_cell{1}.burst_dat{2,1}.start/60 proto_cell{1}.burst_dat{2,1}.stop/60];
c = [ones(1,length(proto_cell{1}.burst_dat{2,1}.start)) zeros(1,length(proto_cell{1}.burst_dat{2,1}.stop))];
[a,b] = sort(burst);
burst = a;
c = c(b);

hold on;
p2 = stairs(burst, c,'linewidth',2);
%++++=====================
% mRNA profile plot
%++++=====================
f = figure;
p1 = stairs(proto_cell{1}.mrna_prof{1}.t/60, proto_cell{1}.mrna_prof{1}.v,'linewidth',2);
hold on;
p2 = stairs(proto_cell{1}.mrna_prof{2,1}.t/60, proto_cell{1}.mrna_prof{2,1}.v,'linewidth',2);

%++++=====================
% Protein profile plot
%++++=====================
f = figure;
p1 = stairs(proto_cell{1}.prot_prof{1}.t/60, proto_cell{1}.prot_prof{1}.v,'linewidth',2);
hold on;
p2 = stairs(proto_cell{1}.prot_prof{2,1}.t/60, proto_cell{1}.prot_prof{2,1}.v,'linewidth',2);


%++++=====================
% Metabolite plot
%++++=====================
f = figure;
p1 = plot(proto_cell{1}.metab_prof.x/60, proto_cell{1}.metab_prof.y(end-p+1:end,:),'linewidth',2);
hold on;
p2 = plot(xlim,[threshold threshold]*2,'--','linewidth',2,'color','k');

% Readjust size 
f.Position = [403 246 560 190];	% Default = [403 246 560 420]: 1) XPos screen, 2) Ypos screen, 3) X size, 4) Y size

% Increase font
f.Children.XAxis.FontSize = 11;
f.Children.YAxis.FontSize = 11;

%++++=====================
plot(proto_cell{1}.mrna_prof{1}.t ,proto_cell{1}.mrna_prof{1}.v)
stairs(proto_cell{1}.mrna_prof{1}.t ,proto_cell{1}.mrna_prof{1}.v,'linewidth',2)



