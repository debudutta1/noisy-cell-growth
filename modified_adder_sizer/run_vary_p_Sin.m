% MONOD'S OBSERVATION - VARY p and Sin
% RUN SIMULATION



p = [1 2 3 5 10];
Sin = 1e3*[5 7.5 11 20 50];
seed_range = 0:4;

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
n = 3;
T = 1.5;
createRec = 0;

same_thresh = 1;
if same_thresh == 1
	threshold = 1.9745e+07;
else
	for i = 1:p
		threshold(i) = 1.9745e+07;
	end
end

%base_dir = 'D:\Debu Simulations\Dec 2019\Data';
base_dir = 'D:\Debu Simulations\Dec 2019\Data\modified_adder_sizer';

%mkdir(strcat(base_dir,'\modified_adder_sizer'));
mkdir(strcat(base_dir,'\var_p_sin'));
%mkdir('D:\Debu Simulations\Dec 2019\Data\var_p_sin\Synced');
mkdir(strcat(base_dir,'\var_p_sin\Synced'));

% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
mkdir(strcat(base_dir,'\var_p_sin\Synced\Sizer'));

for j = 1:length(p)
	for i = 1:length(Sin)
		for k = seed_range
			rng(k);
			y = parallel_growth_sim_v2(gens, Sin(i), n,p(j), T, config, sec_or_min, createRec, threshold, same_thresh);
			%save(strcat('D:\Debu Simulations\Dec 2019\Data\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'.mat'),'y');
			div_durs = y.div_durs;
			div_durs_exp = y.div_durs_exp;
			proto_cell = y.proto_cell;
			cell_rec = y.cell_rec;
			size_bir_exp = y.size_bir_exp;
			size_div_exp = y.size_div_exp;
			sim_vars = struct('Sin',Sin(i),'p',p(j));
			save(strcat(base_dir,'\var_p_sin\Synced\Sizer\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp');
		end
	end
end


% ADDER with no delay - Synced
config = struct('adder',3,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent




mkdir(strcat(base_dir,'\var_p_sin\Synced\Adder'));

for j = 1:length(p)
	for i = 1:length(Sin)
		for k = seed_range
			rng(k);
			y = parallel_growth_sim_v2(gens, Sin(i), n,p(j), T, config, sec_or_min, createRec, threshold, same_thresh);
			%save(strcat('D:\Debu Simulations\Dec 2019\Data\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'.mat'),'y');
			div_durs = y.div_durs;
			div_durs_exp = y.div_durs_exp;
			proto_cell = y.proto_cell;
			cell_rec = y.cell_rec;
			size_bir_exp = y.size_bir_exp;
			size_div_exp = y.size_div_exp;
			sim_vars = struct('Sin',Sin(i),'p',p(j));
			save(strcat(base_dir,'\var_p_sin\Synced\Adder\','var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k),'.mat'),'div_durs','div_durs_exp','proto_cell','sim_vars','config','size_bir_exp','size_div_exp');
			%save(strcat('D:\Debu Simulations\Dec 2019\Data\var_p_sin\Synced\Adder\rng',num2str(k),'\var_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'.mat'),'div_durs','div_durs_exp','proto_cell','cell_rec','sim_vars');
		end
	end
end

