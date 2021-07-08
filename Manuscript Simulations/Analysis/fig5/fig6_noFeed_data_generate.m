Sin = 1e3*[7.5 10 13 20];
p = [3];
seed_range = 0:2;

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;

createRec = 2;

same_thresh = 1;
if same_thresh == 1
	threshold = 2e+07;
else
	for i = 1:p
		threshold(i) = 2e+07;
	end
end

%base_dir = 'D:\Debu Simulations\Dec 2019\Data';
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';

%mkdir(strcat(base_dir,'\modified_adder_sizer'));
mkdir(strcat(base_dir,'\var_p_sin'));
%mkdir('D:\Debu Simulations\Dec 2019\Data\var_p_sin\Synced');
mkdir(strcat(base_dir,'\var_p_sin\fig6_cf_sizer'));

% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent

for k = seed_range
	for j = 1:length(p)
		for i = length(Sin):-1:1
			
			rng(k);
			y = parallel_growth_sim_v2_1(gens, Sin(i), n,p(j), T, config, sec_or_min, createRec, threshold, same_thresh);
			
			div_durs = y.div_durs;
			div_durs_exp = y.div_durs_exp;

			x_dat(k+1,:,i,j) = div_durs_exp;
			
			sim_vars = struct('Sin',Sin(i),'p',p(j),'seed',k);
			
			save(strcat(base_dir,'\var_p_sin\fig6_cf_sizer\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','sim_vars','config');
			
		end
	end
end

%% calculate growth_rate

addpath('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\')
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

save('fig6_cf_sizer_nofeed_data.mat','growth_rate','x_dat','Sin');