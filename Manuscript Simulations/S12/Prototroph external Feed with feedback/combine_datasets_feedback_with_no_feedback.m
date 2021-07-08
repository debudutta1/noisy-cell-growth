% combine data sets of ext Feed with and without feedback

path_loc = 'D:\Debu Simulations\Sep 2020\var_extFeed\';

ld_dat = load(strcat(path_loc,'var_extFeed_fig4_seed_0-6.mat'),'growth_rate','div_durs_compiled');%,'Sin','p','seed_range');
growth_rate_o = ld_dat.growth_rate;
div_durs_compiled_o = ld_dat.div_durs_compiled;
clear('ld_dat');

% =====================
% rename data variables

growth_rate_n = growth_rate;
div_durs_compiled_n = div_durs_compiled;

clear('growth_rate','div_durs_compiled')

% =====================
base_dir = 'D:\Debu Simulations\Sep 2020\var_extFeed_Rev1\';

k3_var_range = [0 1 2.5 5 10 25 50];	% Load 0 data from normal simulation
n_ext = 1;

for k3_var = 1:length(k3_var_range)
	for j = 1:length(feed_rate)
		for n_ext = n_ext_range
			% If no feed, n_ext doesn't matter. Skip redundant simulations
			if feed_rate(j) == 0 & n_ext > 1
				break;
			end
			for i = 1:length(Sin)
				for k = seed_range+1
					
					if k3_var == 1
						div_durs_compiled(:,k,j,n_ext,i,k3_var) = div_durs_compiled_o(:,k,j,n_ext,i);
						growth_rate(k,j,n_ext,i,k3_var) = growth_rate_o(k,j,n_ext,i);
					else
					
						div_durs_compiled(:,k,j,n_ext,i,k3_var) = div_durs_compiled_n(:,k,j,n_ext,i,k3_var-1); 
						growth_rate(k,j,n_ext,i,k3_var) = growth_rate_n(k,j,n_ext,i,k3_var-1);
					end	
				end
			end
		end
	end
end

save(strcat(base_dir,'var_extFeed_rev1_var_k3_withNoFeedback.mat'),'div_durs_compiled','growth_rate','n_ext','k3_var_range');