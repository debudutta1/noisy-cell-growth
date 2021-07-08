%% ANALYSIS ALT_auxFeedSec_noisyTransport

sec_or_min = 0; % sec
n = 3;
p = 3;
n_ext = 0;
seed_range = 0:3;

Sin = 1e3*[1.5 3 4 5 6 6.5 7 9 13];

OVprod_mode = 2; % Decide the scheme of production of secreted metabolite: 1) Only Flux Doubling, 2) Only enzyme Overexpression, 3) Both 

% Noisy Feed rates are fractions to multiply with k3.
feed_rate = [0.5 0.8 0.9 1 1.1 1.2];

secRatio = 0;	% Only secretions

gens = 13; % 2^14 cells = 16384

base_dir = 'D:\Debu Simulations\Sep 2020\var_auxFeed_noisyTransport\';

for k = seed_range+1	%seed_range
	for i = 1:length(Sin)
		for j = length(feed_rate):-1:1
			for l = 1:length(secRatio)
			
				load(strcat(base_dir,'ALT_aux1_Sin',num2str(Sin(i)),'_secR',num2str(secRatio(l)),'_feed',num2str(feed_rate(j)),'_rng',num2str(k-1),'.mat'),'div_durs_exp','size_bir_exp','size_div_exp'); %,'secretion_prof'
				
				div_durs_compiled(:,k,j,i) = div_durs_exp;
				size_bir_compiled(:,:,k,j,i) = size_bir_exp;
				size_div_compiled(:,:,k,j,i) = size_div_exp;
				
			end
		end
	end
end

	parpool(length(seed_range))
	growth_rate = nan(length(seed_range),length(feed_rate),length(Sin));
	reps = 100;
	for j = 1:length(feed_rate)
		for i = 1:length(Sin)
			parfor k = seed_range+1
				tic; growth_rate(k,j,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,i)); toc;
			end
		end
	end



save(strcat(base_dir,'rev2_auxFeed_noisy_dat.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate'); % 'secretion_prof'

%%========================
%% load var_sin_p data
%% load original prototroph var_sin_p data
path_loc = 'C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Monod\';
ld_dat = load(strcat(path_loc,'divD_growR_vary_sin_p_sizer.mat'),'growth_rate');%,'Sin','p','seed_range');
growth_rate_o = ld_dat.growth_rate;
clear('ld_dat');
load(strcat(path_loc,'divD_growR_vary_sin_p_sizer.mat'),'x_dat');%,'Sin','p','seed_range');
%%========================
% Line COLORS
lin_colors = lines(14);

sel = [2 3 4 6 7 8 10 13];	% %% Select the data points for Sin - growth_rate_o
sel_dat = 2:length(Sin);

plot_indices = [1:4 6];	% feed_ratio

% growth_rate(k,j,i) => k-Seed, j - feed_rate, i - Sin
fg = figure; hold on;
% Plot in black dotted growth rates of p = 2, and p = 3 prototroph	
j = 3; %=>p=3;
plot(Sin(sel_dat), reshape(mean(growth_rate_o(:,sel,j),1),length(Sin(sel_dat)),[]),'--o','Displayname',strjoin({'p =', num2str(j),'\newline'},' '),'linewidth',2,'color',[.4 .4 .4])
j = 2; %=>p=2;
plot(Sin(sel_dat), reshape(mean(growth_rate_o(:,sel,j),1),length(Sin(sel_dat)),[]),'--o','Displayname',strjoin({'p =', num2str(j)},' '),'linewidth',2,'color','k')

%for j = 2:length(feed_rate)
for j = 1:length(plot_indices)
	plot(Sin(sel_dat), reshape(mean(growth_rate(:,plot_indices(j),sel_dat),1),length(Sin(sel_dat)),[]),'-o','Displayname',strjoin({'Feed ratio =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',2,'color',lin_colors(j,:))
end
%title({'Variation of growth rate of auxotroph with varying feed','and varying Substrate flux',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux (per sec)')
%legend('location','east')%,'NumColumns',3);
%legend('location','southoutside','NumColumns',3);
legend('location','east','NumColumns',1);

xlim([2700 13300])
ylim([1 1.37])
%xlim([1250 13250])

%addpath('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\MATLAB file exchange\Break Axis\breakyaxis\')
%breakyaxis([.33 .58]);

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = 'fig6a_growthR_Sin_Noisyfeed';
%savefig(strcat(name,'.fig'));
%print(name,'-dpng','-r300');
set(fg, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r300', '-transparent', '-painters')


%%========================
%%========================

sec_or_min = 0; % sec
n = 3;
p = 3;
n_ext = 0;

% define seed_range before running script
if exist('seed_range') == 0
	seed_range = 0:3;	
end

%Sin = 1e3*[1.5 4 5.5 6.5 7 7.5 10 20];
%Sin = 1e3*[1.5 3 4 5 6 6.5 7 9 13];
Sin = 1e3*[4 5 6 6.5 7 9];

% Noisy Feed rates are fractions to multiply with k3.
feed_rate = [0.5 0.8 0.9 1 1.1 1.2 2];

secRatio = 0.5;	% Only 0.5 to compare with the non-noisy CF case

OVprod_mode = 1; % Decide the scheme of production of secreted metabolite: 1) Only Flux Doubling, 2) Only enzyme Overexpression, 3) Both 

% Exponential distribution D1 to derive the duration of transcription ON time
t_ON = 4; % On Average transcription ON time in minutes
% Exponential distribution D2 to derive the duration of transcription OFF time
t_OFF = 2.4; % On Average transcription OFF time in minutes

genExp.t_ON = t_ON;
genExp.t_OFF = t_OFF;
% NO Overexpression
genExp.t_ON_ov = t_ON;
genExp.t_OFF_ov = t_OFF;

threshold = 1e7;
same_thresh = 1;

%addpath('D:\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\');

% NEW ADDER with no delay - Synced
%config = struct('adder',3,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
% NEW SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent

gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
%createRec = 0;
createRec = 2;	% Detailed data not output, but computed to extract important parameters

% ALT version: Secretion from final metabolite 
cell_type1 = [1, 3, 2, 3];	%[1, 2, 2, 3];	%[type, sec_n, sec_p, aux_p]

base_dir = 'D:\Debu Simulations\Sep 2020\';
base_dir = strcat(base_dir,'var_auxFeedSecrete_noisyTransport\');

for k = seed_range+1	%seed_range
	for i = 1:length(Sin)
		for j = length(feed_rate):-1:1
			for l = 1:length(secRatio)
			
				load(strcat(base_dir,'ALT_aux1_Sin',num2str(Sin(i)),'_secR',num2str(secRatio(l)),'_feed',num2str(feed_rate(j)),'_rng',num2str(k-1),'.mat'),'div_durs_exp','size_bir_exp','size_div_exp'); %,'secretion_prof'
				
				div_durs_compiled(:,k,j,i) = div_durs_exp;
				size_bir_compiled(:,:,k,j,i) = size_bir_exp;
				size_div_compiled(:,:,k,j,i) = size_div_exp;
				
			end
		end
	end
end

	parpool(length(seed_range))
	growth_rate = nan(length(seed_range),length(feed_rate),length(Sin));
	reps = 100;
	for j = 1:length(feed_rate)
		for i = 1:length(Sin)
			parfor k = seed_range+1
				tic; growth_rate(k,j,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,i)); toc;
			end
		end
	end

save(strcat(base_dir,'rev2_auxFeedSec_noisy_dat.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate'); % 'secretion_prof'


%%================================
nofeed = load('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Monod\divD_growR_vary_sin_p_Sizer.mat');

% ExtFeed Prototroph data with Feed
extfeed = load(strcat('D:\Debu Simulations\Sep 2020\var_extFeed\var_extFeed_fig4.mat'));

% AuxFeed
auxfeed = load('D:\Debu Simulations\Sep 2020\var_auxFeed\auxFeed_fig5_dat2.mat','div_durs_compiled','growth_rate'); 
auxfeed.Sin = 1e3*[1.5 3 4 5 6 6.5 7 9 13];
%%================================
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
plot_indices = [1:4 6];	% feed_rate
%for j = 2:length(feed_rate)
for j = 1:length(plot_indices)
	plot(Sin(sel_dat), reshape(mean(growth_rate(:,plot_indices(j),sel_dat),1),length(Sin(sel_dat)),[]),'-o','Displayname',strjoin({'Feed =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',2,'color',lin_colors(j,:))
end
%title({'Variation of growth rate of auxotroph with varying feed','and varying Substrate flux',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux (per sec)')
%legend('location','east')%,'NumColumns',3);
%legend('location','southoutside','NumColumns',3);
legend('location','east','NumColumns',2);

f.Position(3:4) = f.Position(3:4)*0.75;
f.Children(2).XLabel.FontSize = 12;
f.Children(2).YLabel.FontSize = 12;

name = 'fig6b_growthR_Sin_Noisyfeed_Sec';
set(f, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r300', '-transparent', '-painters')

%%================================


% Line COLORS
lin_colors = lines(14);

sel = [2 3 4 6 7 8 10 13];	% %% Select the data points for Sin - growth_rate_o
sel_dat = 2:length(Sin);

plot_indices = [1:4 6];	% feed_ratio

% growth_rate(k,j,i) => k-Seed, j - feed_rate, i - Sin
fg = figure; hold on;
% Plot in black dotted growth rates of p = 2, and p = 3 prototroph	
j = 3; %=>p=3;
plot(Sin(sel_dat), reshape(mean(growth_rate_o(:,sel,j),1),length(Sin(sel_dat)),[]),'--o','Displayname',strjoin({'p =', num2str(j),'\newline'},' '),'linewidth',2,'color',[.4 .4 .4])
j = 2; %=>p=2;
plot(Sin(sel_dat), reshape(mean(growth_rate_o(:,sel,j),1),length(Sin(sel_dat)),[]),'--o','Displayname',strjoin({'p =', num2str(j)},' '),'linewidth',2,'color','k')

%for j = 2:length(feed_rate)
for j = 1:length(plot_indices)
	plot(Sin(sel_dat), reshape(mean(growth_rate(:,plot_indices(j),sel_dat),1),length(Sin(sel_dat)),[]),'-o','Displayname',strjoin({'Feed ratio =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',2,'color',lin_colors(j,:))
end
%title({'Variation of growth rate of auxotroph with varying feed','and varying Substrate flux',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux (per sec)')
%legend('location','east')%,'NumColumns',3);
%legend('location','southoutside','NumColumns',3);
legend('location','east','NumColumns',1);

xlim([2700 13300])
ylim([1 1.37])
%xlim([1250 13250])

%addpath('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\MATLAB file exchange\Break Axis\breakyaxis\')
%breakyaxis([.33 .58]);

fg.Position(3:4) = fg.Position(3:4)*0.75;
fg.Children(2).XLabel.FontSize = 12;
fg.Children(2).YLabel.FontSize = 12;

name = 'fig6a_growthR_Sin_Noisyfeed';
%savefig(strcat(name,'.fig'));
%print(name,'-dpng','-r300');
set(fg, 'Color', 'none');
export_fig(name, '-pdf', '-png', '-r300', '-transparent', '-painters')
