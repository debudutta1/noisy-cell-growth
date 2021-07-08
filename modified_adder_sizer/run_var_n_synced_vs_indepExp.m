%% Simulation of varying n and Independent expression vs Synced

p = [1, 2, 3];
Sin = 1e3*[5 7.5 11 20 50];

seed_range = 0:2;

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = [1 2 3];
T = 1.5;

%createRec = 0;
createRec = 2;	% Detailed data not output, but computed to extract important parameters

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
mkdir(strcat(base_dir,'\var_n_sync_indep'));
mkdir(strcat(base_dir,'\var_n_sync_indep\Synced'));

% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
mkdir(strcat(base_dir,'\var_n_sync_indep\Synced\Sizer'));

for k = seed_range
	for j = 1:length(p)
		for l = 1:length(n)
			for i = length(Sin):-1:1
				
				rng(k);
				y = parallel_growth_sim_v2_1(gens, Sin(i), n(l),p(j), T, config, sec_or_min, createRec, threshold, same_thresh);
				
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
					
				sim_vars = struct('Sin',Sin(i),'p',p(j),'seed',k,'n',n(l));
				
				save(strcat(base_dir,'\var_n_sync_indep\Synced\Sizer\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'_n',num2str(n(l)),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced');
			end
		end
	end
end


%mkdir(strcat(base_dir,'\modified_adder_sizer'));
mkdir(strcat(base_dir,'\var_n_sync_indep'));
mkdir(strcat(base_dir,'\var_n_sync_indep\Independent'));

% SIZER with Independent Expression
config = struct('adder',2,'reset_S',0,'singOper',0);	% singOper 2 synced, 1 delayed, 0 independent
mkdir(strcat(base_dir,'\var_n_sync_indep\Independent\Sizer'));

for k = seed_range
	for j = 1:length(p)
		for l = 1:length(n)
			for i = length(Sin):-1:1
				
				rng(k);
				y = parallel_growth_sim_v2_1(gens, Sin(i), n(l),p(j), T, config, sec_or_min, createRec, threshold, same_thresh);
				
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
					
				sim_vars = struct('Sin',Sin(i),'p',p(j),'seed',k,'n',n(l));
				
				save(strcat(base_dir,'\var_n_sync_indep\Independent\Sizer\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'_n',num2str(n(l)),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced');
			end
		end
	end
end