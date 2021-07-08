%% Analysis SECRETION PROTOTROPH

addpath('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned')

%Sin = 1e3*[1.5 4 5.5 6.5 7 7.5 10 20];
Sin = 1e3*[1.5 5 6 6.5 7 9 13];
p = [3 1 5];
%seed_range = 0:1;
seed_range = 0:2;

secRatio = [0 0.05 0.1 0.25 0.5 0.75];

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;

%createRec = 0;
createRec = 2;	% Detailed data not output, but computed to extract important parameters

same_thresh = 1;
if same_thresh == 1
	threshold = 1e+07;
else
	for i = 1:p
		threshold(i) = 1e+07;
	end
end

base_dir = 'D:\Debu Simulations\Sep 2020\var_protSec\';

% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent

size_bir_compiled = nan(2^gens-1,length(p),length(seed_range),length(p),length(secRatio),length(Sin));
size_div_compiled = nan(2^gens-1,length(p),length(seed_range),length(p),length(secRatio),length(Sin));

for k = seed_range+1
	for j = 1:length(p)
		for i = length(Sin):-1:1
			for l = 1:length(secRatio)
				
				load(strcat(base_dir,'var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k-1),'_secRatio',num2str(secRatio(l)),'.mat'),'div_durs_exp','sim_vars','size_bir_exp','size_div_exp');
				
				div_durs_compiled(:,k,j,l,i) = div_durs_exp;
%				size_bir_compiled(:,:,k,j,l,i) = size_bir_exp;
%				size_div_compiled(:,:,k,j,l,i) = size_div_exp;

			end
		end
	end
end

% Obtain Growth rate of cells from dataset
parpool(seed_range(end)+1)
growth_rate = nan(length(seed_range),length(p),length(secRatio),length(Sin));
reps = 100;
for j = 1:length(p)
	for i = 1:length(Sin)
		for l = 1:length(secRatio)
			parfor k = seed_range+1
				tic; growth_rate(k,j,l,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,l,i)); toc;
			end
		end
	end
end

save(strcat(base_dir,'var_protSec_fig4.mat'),'div_durs_compiled','growth_rate');
%save(strcat(base_dir,'var_protSec_fig4.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate');
