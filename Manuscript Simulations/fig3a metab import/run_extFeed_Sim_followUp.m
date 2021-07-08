% Run oneFeed Simulations

sec_or_min = 0; % sec
n = 3;
p = 3; % [2,3,5];
n_ext_range = [1:p];
seed_range = 0:4;
%feed_rate = [0:2500:10000];
%feed_rate = [12500 15000];

%feed_rate = [0:2500:15000];
feed_rate = [625 1250];

t_ON = 6; 
%t_OFF = 9;
t_OFF = 37;
threshold = 2e7;
same_thresh = 1;

addpath('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\');

% NEW ADDER with no delay - Synced
%config = struct('adder',3,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
% NEW SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent

%Sin = 1e3*[5 7.5 11 20 50];
%Sin = 1e3*[7.5 11 13 20];
Sin = 1e3*[7.5 10 13 20 50];

gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
createRec = 0;

base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';
%mkdir(strcat(base_dir,'\var_extFeed'));
%mkdir(strcat(base_dir,'\var_extFeed\sizer'));
%base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_extFeed\sizer\';
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer\var_extFeed\sizer_toff_37\';

for k = seed_range
	for j = 1:length(feed_rate)
		for n_ext = n_ext_range
			% If no feed, n_ext doesn't matter. Skip redundant simulations
			if feed_rate(j) == 0 & n_ext > 1
				break;
			end
			for i = 1:length(Sin)
				
				rng(k)
				display(strcat('seed_',num2str(k),'feed_',num2str(j),'n_ext_',num2str(n_ext),'Sin_',num2str(i)))
				if isfile(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat'))
					display(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat'))
					continue;
				end
				
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
				
				save(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','sim_vars','config','size_bir','size_div','size_bir_exp','size_div_exp');
			end
		end
	end
end