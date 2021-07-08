% MONOD'S OBSERVATION - VARY p and Sin
% RUN SIMULATION

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

base_dir = 'D:\Debu Simulations\Sep 2020\';
mkdir(strcat(base_dir,'\var_protSec'));
base_dir = 'D:\Debu Simulations\Sep 2020\var_protSec\';

% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent

for k = seed_range
	for j = 1:length(p)
		for i = length(Sin):-1:1
			for l = 1:length(secRatio)
				
				if exist(strcat(base_dir,'var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'_secRatio',num2str(secRatio(l)),'.mat')) == 2
					continue
				end		
				
				rng(k);
				y = parallel_growth_sim_protSec(gens, Sin(i), n,p(j), T, config, sec_or_min, createRec, threshold, same_thresh, secRatio(l));
				
				div_durs = y.div_durs;
				div_durs_exp = y.div_durs_exp;
				proto_cell = y.proto_cell;
				cell_rec = y.cell_rec;
				size_bir = y.size_bir;
				size_div = y.size_div;
				size_bir_exp = y.size_bir_exp;
				size_div_exp = y.size_div_exp;
					
				mrna_produced = y.mrna_produced;
				prot_produced = y.prot_produced;
				mrna_beg = y.mrna_beg;
				mrna_end = y.mrna_end;
				prot_beg = y.prot_beg;
				prot_end = y.prot_end;
					
				sim_vars = struct('Sin',Sin(i),'p',p(j),'seed',k,'secRatio',secRatio(l));
				
				save(strcat(base_dir,'var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'_secRatio',num2str(secRatio(l)),'.mat'),'div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp');
				%save(strcat(base_dir,'var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'_secRatio',num2str(secRatio(l)),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced');
			end
		end
	end
end