% Gene Expression parameter var and thereshold var

Sin = 10000e3;
%p = 1;
%p = 3;
p = 5;

t_off_range = 37;
t_on_range = 6;
threshold_range = [1e4, 5e4, 1e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9];

seed_range = 0:2;

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;
createRec = 0;

same_thresh = 1;

base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';
mkdir(strcat(base_dir,'\var_threshold'));
mkdir(strcat(base_dir,'\var_threshold\sizer'));

% NEW SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent


for i = 1:length(threshold_range)
	for j = 1:length(t_off_range)
		for l = 1:length(t_on_range)		
			for k = seed_range
				rng(k);
				y = parallel_growth_sim_v3(gens, Sin,n,p, T, config, sec_or_min, createRec, threshold_range(i), same_thresh, t_on_range(l), t_off_range(j));
				%save(strcat('D:\Debu Simulations\Dec 2019\Data\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'.mat'),'y');
				div_durs = y.div_durs;
				div_durs_exp = y.div_durs_exp;
				proto_cell = y.proto_cell;
				cell_rec = y.cell_rec;
				size_bir = y.size_bir;
				size_div = y.size_div;
				size_bir_exp = y.size_bir_exp;
				size_div_exp = y.size_div_exp;
				sim_vars = struct('threshold',threshold_range(i),'t_on',t_on_range(l), 't_off',t_off_range(j));
				save(strcat(base_dir,'\var_threshold\sizer\','var_thr',num2str(threshold_range(i)),'_ton',num2str(t_on_range(l)),'_toff',num2str(t_off_range(j)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir','size_div','size_bir_exp','size_div_exp');
			end
		end
	end
end


% NEW ADDER with no delay - Synced
config = struct('adder',3,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
mkdir(strcat(base_dir,'\var_threshold\adder'));

for i = 1:length(threshold_range)
	for j = 1:length(t_off_range)
		for l = 1:length(t_on_range)		
			for k = seed_range
				rng(k);
				y = parallel_growth_sim_v3(gens, Sin,n,p, T, config, sec_or_min, createRec, threshold_range(i), same_thresh, t_on_range(l), t_off_range(j));
				%save(strcat('D:\Debu Simulations\Dec 2019\Data\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'.mat'),'y');
				div_durs = y.div_durs;
				div_durs_exp = y.div_durs_exp;
				proto_cell = y.proto_cell;
				cell_rec = y.cell_rec;
				size_bir = y.size_bir;
				size_div = y.size_div;
				size_bir_exp = y.size_bir_exp;
				size_div_exp = y.size_div_exp;
				sim_vars = struct('threshold',threshold_range(i),'t_on',t_on_range(l), 't_off',t_off_range(j));
				save(strcat(base_dir,'\var_threshold\adder\','var_thr',num2str(threshold_range(i)),'_ton',num2str(t_on_range(l)),'_toff',num2str(t_off_range(j)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir','size_div','size_bir_exp','size_div_exp');
			end
		end
	end
end