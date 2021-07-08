% Analysis Prototroph + Secrete

addpath('C:\Users\sysbio_admin\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer')

Sin = 1e3*[7.5 10 13 20 50];
%Sin = 1e3*[5 7.5 11 20 50];
% p = [1 3 5];
p = 3;
seed_range = 0:2;

secRatio = [0 0.05 0.1 0.25 0.5 0.75];
%secRatio = [0.05 0.1 0.25 0.5 0.75];

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;

base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';

size_bir_compiled = nan(2^gens-1,p,length(seed_range),length(p),length(secRatio),length(Sin));
size_div_compiled = nan(2^gens-1,p,length(seed_range),length(p),length(secRatio),length(Sin));

for k = seed_range+1
	for j = 1:length(p)
		for i = 1:length(Sin)
			for l = 1:length(secRatio)
				load(strcat(base_dir,'\protSec\Synced\Sizer\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k-1),'_secRatio',num2str(secRatio(l)),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced');
				
				div_durs_compiled(:,k,j,l,i) = div_durs_exp;
				size_bir_compiled(:,:,k,j,l,i) = size_bir_exp;
				size_div_compiled(:,:,k,j,l,i) = size_div_exp;
			end
		end
	end
end

%% Obtain Growth rate of cells from dataset
parpool(seed_range(end)+1)
growth_rate = nan(length(seed_range),length(p),length(secRatio),length(Sin));
reps = 100;
for i = 1:length(Sin)
	for j = 1:length(p)
		for l = 1:length(secRatio)
			parfor k = seed_range+1
				tic; growth_rate(k,j,l,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,l,i)); toc;
			end
		end
	end
end
%save(strcat(base_dir,'var_extFeed_growR.mat'),'growth_rate');

save(strcat(base_dir,'\protSec\Synced\Sizer\','var_protSec_fig4_p3.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate');