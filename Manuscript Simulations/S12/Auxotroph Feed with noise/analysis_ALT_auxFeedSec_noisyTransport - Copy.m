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


%% ANALYSIS ALT_auxFeedSec

sec_or_min = 0; % sec
n = 3;
p = 3;
n_ext = 0;

seed_range = 0:4;
Sin = 1e3*[7.5 10 13 20];
feed_rate = [7500:2500:15000];
secRatio = [0.3:0.1:0.7];
gens = 13; % 2^14 cells = 16384

addpath('C:\Users\Dibyendu_Lab\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\');


base_dir = 'E:\Debu_Work\Dec 2019\var_ALT_auxFeedSec\sizer\1_fluxDouble\';
load(strcat(base_dir,'1_fluxDouble_ALT_Dat.mat'));
load(strcat(base_dir,'auxFeedSecrete_growR.mat'),'growth_rate');


base_dir = 'E:\Debu_Work\Dec 2019\var_ALT_auxFeedSec\sizer\2_enzDouble_par6\';
load(strcat(base_dir,'2_enzDouble_par6_ALT_Dat.mat'));
load(strcat(base_dir,'auxFeedSecrete_growR.mat'),'growth_rate');


base_dir = 'E:\Debu_Work\Dec 2019\var_ALT_auxFeedSec\sizer\2_enzDouble_par7\';
load(strcat(base_dir,'2_enzDouble_par7_ALT_Dat.mat'));
load(strcat(base_dir,'auxFeedSecrete_growR.mat'),'growth_rate');


base_dir = 'E:\Debu_Work\Dec 2019\var_ALT_auxFeedSec\sizer\2_enzDouble_par2\';
load(strcat(base_dir,'2_enzDouble_par2_ALT_Dat.mat'));
load(strcat(base_dir,'auxFeedSecrete_growR.mat'),'growth_rate');


