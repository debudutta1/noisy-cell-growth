%% ANALYSIS ALT_auxFeedSec

sec_or_min = 0; % sec
n = 3;
p = 3;
n_ext = 0;
seed_range = 0:3;

%Sin = 1e3*[1.5 5 6 6.5 7 9 13];
Sin = 1e3*[4 5 6 6.5 7 9];

OVprod_mode = 1; % Decide the scheme of production of secreted metabolite: 1) Only Flux Doubling, 2) Only enzyme Overexpression, 3) Both 

%feed_rate = [312 625 1250 2500 5000 7500];% 12500];
feed_rate = 1e3*[4 5 5.5 6 6.5 7];% 12500];

%secRatio = 0.5;
secRatio = [0.3:0.1:0.7];

gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
%createRec = 0;
createRec = 2;	% Detailed data not output, but computed to extract important parameters

% ALT version: Secretion from final metabolite 
cell_type1 = [1, 3, 2, 3];	%[1, 2, 2, 3];	%[type, sec_n, sec_p, aux_p]

addpath('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned');

base_dir = 'D:\Debu Simulations\Sep 2020\';
if OVprod_mode == 1
%	mkdir(strcat(base_dir,'var_ALT_auxFeedSec\1_fluxDouble'));
	base_dir = strcat(base_dir,'var_ALT_auxFeedSec\1_fluxDouble\');
elseif OVprod_mode == 2
%	mkdir(strcat(base_dir,'var_ALT_auxFeedSec\2_enzDouble_par2'));
	base_dir = strcat(base_dir,'var_ALT_auxFeedSec\2_enzDouble_par2\');
elseif OVprod_mode == 3
%	mkdir(strcat(base_dir,'var_ALT_auxFeedSec\3_FluxEnzDoub'));
	base_dir = strcat(base_dir,'var_ALT_auxFeedSec\3_FluxEnzDoub\');
end

for k = seed_range+1	%seed_range
	for i = 1:length(Sin)
		for j = length(feed_rate):-1:1
			for l = 1:length(secRatio)
				
				load(strcat(base_dir,'ALT_aux1_Sin',num2str(Sin(i)),'_secR',num2str(secRatio(l)),'_feed',num2str(feed_rate(j)),'_rng',num2str(k-1),'.mat'),'div_durs_exp','size_bir_exp','size_div_exp','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced','sec_end'); %,'secretion_prof'
				
				div_durs_compiled(:,k,j,l,i) = div_durs_exp;
				size_bir_compiled(:,:,k,j,l,i) = size_bir_exp;
				size_div_compiled(:,:,k,j,l,i) = size_div_exp;
				
%				sec_prof_compiled(:,k,j,l,i) = secretion_prof;
				sec_end_compiled(:,k,j,l,i) = sec_end;
				
				mrna_prod_compiled(:,:,:,k,j,l,i) = mrna_produced;
				prot_prod_compiled(:,:,:,k,j,l,i) = prot_produced;
				
				mrna_beg_compiled(:,:,:,k,j,l,i) = mrna_beg;
				mrna_end_compiled(:,:,:,k,j,l,i) = mrna_end;
				prot_beg_compiled(:,:,:,k,j,l,i) = prot_beg;
				prot_end_compiled(:,:,:,k,j,l,i) = prot_end;
				
			end
		end
	end
end

	parpool(length(seed_range))
	growth_rate = nan(length(seed_range),length(feed_rate),length(secRatio),length(Sin));
	reps = 100;
	for l = 1:length(secRatio)
		for j = 1:length(feed_rate)
			for i = 1:length(Sin)
				parfor k = seed_range+1
					tic; growth_rate(k,j,l,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,l,i)); toc;
				end
			end
		end
	end



folder_n = '1_fluxDouble_ALT';
%folder_n = '2_enzDouble_par6_ALT';
%folder_n = '2_enzDouble_par7_ALT';
%save(strcat(base_dir,folder_n,'_Dat.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','mrna_beg_compiled','mrna_end_compiled','mrna_prod_compiled','prot_beg_compiled','prot_end_compiled','prot_prod_compiled','sec_end_compiled','growth_rate'); % 'secretion_prof'
save(strcat(base_dir,folder_n,'_Dat_v2.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','mrna_beg_compiled','mrna_end_compiled','mrna_prod_compiled','prot_beg_compiled','prot_end_compiled','prot_prod_compiled','sec_end_compiled','growth_rate','Sin','feed_rate','secRatio'); % 'secretion_prof'


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


