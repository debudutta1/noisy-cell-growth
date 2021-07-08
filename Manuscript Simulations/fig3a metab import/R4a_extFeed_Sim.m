% Run oneFeed Simulations

sec_or_min = 0; % sec
n = 3;
p = 3; % [2,3,5];
%n_ext_range = [1:p];
n_ext_range = [1 p];
feed_rate = [0 312 625 1250:1250:5000 7500 12500];
%feed_rate = [0:2500:10000];

%seed_range = 0:2;
seed_range = 3:4;

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

%Sin = 1e3*[5 7.5 11 20 50];
%Sin = 1e3*[3 5 7 9 11 13 20 50];
%Sin = 1e3*[1.5 4 5.5 6.5 7 7.5 10 20];
Sin = 1e3*[1.5 5 6 6.5 7 9 13];

gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
createRec = 0;

base_dir = 'D:\Debu Simulations\Sep 2020\';
mkdir(strcat(base_dir,'\var_extFeed'));
base_dir = 'D:\Debu Simulations\Sep 2020\var_extFeed\';

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