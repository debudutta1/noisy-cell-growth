% MONOD'S OBSERVATION - VARY p and Sin
% RUN SIMULATION



p = [1, 2, 3, 4, 5, 10, 11, 14, 15, 19, 20, 24, 25];
Sin = 1e3*[3 6 7 13];

seed_range = 0:2;

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;	% initial T in hrs

createRec = 0;
%createRec = 2;	% Detailed data not output, but computed to extract important parameters

same_thresh = 1;
if same_thresh == 1
	threshold = 1e+07;
else
	for i = 1:p
		threshold(i) = 1e+07;
	end
end

base_dir = 'D:\Debu Simulations\Sep 2020\';
mkdir(strcat(base_dir,'\fig5d_var_p_sin'));
base_dir = 'D:\Debu Simulations\Sep 2020\fig5d_var_p_sin\';


% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
fold_name = 'Sizer';
mkdir(strcat(base_dir,fold_name));


for k = seed_range
	for j = 1:length(p)
		for i = length(Sin):-1:1
			
			if exist(strcat(base_dir,fold_name,'\var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat')) == 2
				continue
			end
			
			rng(k);
			y = parallel_growth_sim_v2_1(gens, Sin(i), n,p(j), T, config, sec_or_min, createRec, threshold, same_thresh);
			
			div_durs = y.div_durs;
			div_durs_exp = y.div_durs_exp;
			proto_cell = y.proto_cell;
			cell_rec = y.cell_rec;
			size_bir = y.size_bir;
			size_div = y.size_div;
			size_bir_exp = y.size_bir_exp;
			size_div_exp = y.size_div_exp;
			
			sim_vars = struct('Sin',Sin(i),'p',p(j),'seed',k);
			
			if createRec == 2
				mrna_produced = y.mrna_produced;
				prot_produced = y.prot_produced;
				mrna_beg = y.mrna_beg;
				mrna_end = y.mrna_end;
				prot_beg = y.prot_beg;
				prot_end = y.prot_end;
				
				save(strcat(base_dir,fold_name,'\var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced');
			else
				save(strcat(base_dir,fold_name,'\var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp');
			end
		end
	end
end


