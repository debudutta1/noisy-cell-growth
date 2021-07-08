% Analysis of Fig 5d var_sin_p simulations

% Parameters
p = [1, 2, 3, 4, 5, 10, 11, 14, 15, 19, 20, 24, 25];
Sin = 1e3*[3 6 7 13];

seed_range = 0:2;

% Compile Data from Raw Simulation data
base_dir = 'D:\Debu Simulations\Sep 2020\fig5d_var_p_sin\';

% SIZER
fold_name = 'Sizer';
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
%parpool('local',3);
growth_rate = nan(length(seed_range),length(Sin),length(p));
reps = 100;
for j = 1:length(p)
	for i = 1:length(Sin)
		parfor k = seed_range+1
			tic; growth_rate(k,i,j) = exp_grow_rate(reps, gens, x_dat(k,:,i,j)); toc;
		end
	end
end

save(strcat('fig5d_dat_',fold_name,'.mat'),'growth_rate', 'x_dat','Sin','p','seed_range');