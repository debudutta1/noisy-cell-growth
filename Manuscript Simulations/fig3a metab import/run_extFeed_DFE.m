% Run oneFeed Simulations

sec_or_min = 0; % sec
n = 3;
p = 3; % [2,3,5];
%n_ext_range = [1:p];
n_ext_range = [1];

feed_rate = linspace(0,8000,20);
%feed_rate = [0 312 625 1250:1250:5000 7500 12500];
%feed_rate = [0:2500:10000];

seed_range = 0:5;
%seed_range = 3:4;

% Exponential distribution D1 to derive the duration of transcription ON time
t_ON = 4; % On Average transcription ON time in minutes
% Exponential distribution D2 to derive the duration of transcription OFF time
t_OFF = 2.4; % On Average transcription OFF time in minutes

threshold = 1e+07;
same_thresh = 1;

addpath('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned')

% NEW ADDER with no delay - Synced
%config = struct('adder',3,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
% NEW SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent

Sin = 1e3*[9];

gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
createRec = 0;

base_dir = 'D:\Debu Simulations\Sep 2020\';
mkdir(strcat(base_dir,'\extFeed_DFE'));
base_dir = 'D:\Debu Simulations\Sep 2020\extFeed_DFE\';

for j = 1:length(feed_rate)
	for n_ext = n_ext_range
		% If no feed, n_ext doesn't matter. Skip redundant simulations
		if feed_rate(j) == 0 & n_ext > 1
			break;
		end
		for i = 1:length(Sin)
			for k = seed_range

				if exist(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat')) == 2
					continue
				end				
					
				rng(k)
				%y = parallel_growth_sim_extFeed(gens, Sin(i), n, p, T, config, sec_or_min, createRec, n_ext);
				y = parallel_growth_sim_extFeed(gens, Sin(i), n, p, T, config, sec_or_min, createRec, threshold, same_thresh, t_ON, t_OFF, n_ext, feed_rate(j));

				div_durs = y.div_durs;
				div_durs_exp = y.div_durs_exp;
				proto_cell = y.proto_cell;
				cell_rec = y.cell_rec;
                size_bir = y.size_bir;
				size_div = y.size_div;
				size_bir_exp = y.size_bir_exp;
				size_div_exp = y.size_div_exp;
				sim_vars = struct('Sin',Sin(i),'p',p,'n_ext',n_ext,'feed_rate',feed_rate(j),'seed',k);
				
				save(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat'),'div_durs_exp','sim_vars','config','size_bir_exp','size_div_exp');
				%save(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','sim_vars','config','size_bir','size_div','size_bir_exp','size_div_exp');
			end
		end
	end
end

%% =======================================

%laptop
addpath('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned')

addpath('C:\Google Drive\MATLAB File Exchange\')
addpath('C:\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\MATLAB file exchange\export_fig-master\')

% ANALYSIS
sec_or_min = 0; % sec
n = 3;
p = 3; % [2,3,5];
n_ext_range = [1];

feed_rate = linspace(0,8000,20);

seed_range = 0:5;

% Exponential distribution D1 to derive the duration of transcription ON time
t_ON = 4; % On Average transcription ON time in minutes
% Exponential distribution D2 to derive the duration of transcription OFF time
t_OFF = 2.4; % On Average transcription OFF time in minutes

threshold = 1e+07;
same_thresh = 1;

addpath('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned')

Sin = 1e3*[9];

gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
createRec = 0;

%base_dir = 'D:\Debu Simulations\Sep 2020\extFeed_DFE\';
base_dir = 'D:\Sim_data\extFeed_DFE\';

for j = 1:length(feed_rate)
	for n_ext = n_ext_range
		% If no feed, n_ext doesn't matter. Skip redundant simulations
		if feed_rate(j) == 0 & n_ext > 1
			break;
		end
		for i = 1:length(Sin)
			for k = seed_range+1
				load(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k-1),'.mat'),'div_durs_exp','sim_vars','config','size_bir_exp','size_div_exp');

				div_durs_compiled(:,k,j,n_ext,i) = div_durs_exp;
				size_bir_compiled(:,:,k,j,n_ext,i) = size_bir_exp;
				size_div_compiled(:,:,k,j,n_ext,i) = size_div_exp;
				
			end
		end
	end
end

% Obtain Growth rate of cells from dataset
%parpool(length(seed_range))
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

%save(strcat(base_dir,'var_extFeed_fig4.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate');
save(strcat(base_dir,'extFeed_DFE.mat'),'div_durs_compiled','growth_rate');
