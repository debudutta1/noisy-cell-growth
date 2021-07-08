% MONOD'S OBSERVATION - VARY p and Sin
% RUN SIMULATION



%p = [1 2 3 5 10];
%p = [ 2 3 5 10];
Sin = 1e3*[7.5 10 13 20];
%seed_range = 0:4;

%p = [1, 2, 3, 5, 10, 11, 14, 15, 19, 20, 24, 25];
%Sin = 1e3*[7.5 11 20 50];

%p = [11, 14, 15, 19, 20, 24, 25];
%Sin = 1e3*[5];
p = [2 3 5];
seed_range = 0:2;

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;

%createRec = 0;
createRec = 2;	% Detailed data not output, but computed to extract important parameters

ratio_p = [0.3 0.25 0.15];	% for 2 3 5
same_thresh = 0;

%base_dir = 'D:\Debu Simulations\Dec 2019\Data';
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';

%mkdir(strcat(base_dir,'\modified_adder_sizer'));
mkdir(strcat(base_dir,'\var_p_sin_mul_thresh_mul_thresh'));
%mkdir('D:\Debu Simulations\Dec 2019\Data\var_p_sin_mul_thresh\Synced');
mkdir(strcat(base_dir,'\var_p_sin_mul_thresh\Synced'));

% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
mkdir(strcat(base_dir,'\var_p_sin_mul_thresh\Synced\Sizer'));

% SIZER with Independent expression
%config = struct('adder',2,'reset_S',0,'singOper',0);	% singOper 2 synced, 1 delayed, 0 independent
%mkdir(strcat(base_dir,'\var_p_sin_mul_thresh\Synced\Sizer'));

for k = seed_range
	for j = 1:length(p)
		if same_thresh == 1
			threshold = 2e+07;
		else
			for i1 = 1:p(j)
				threshold(i1) = 2e+07 * (1 - ratio_p(j)*(i1-1));
			end
		end
		
		for i = length(Sin):-1:1
						
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
				
			mrna_produced = y.mrna_produced;
			prot_produced = y.prot_produced;
			mrna_beg = y.mrna_beg;
			mrna_end = y.mrna_end;
			prot_beg = y.prot_beg;
			prot_end = y.prot_end;
				
			sim_vars = struct('Sin',Sin(i),'p',p(j),'seed',k,'threshold',threshold);
			
			save(strcat(base_dir,'\var_p_sin_mul_thresh\Synced\Sizer\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced');
		end
	end
end