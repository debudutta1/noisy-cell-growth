function bac = cell_growth_fn_v2(n,p,T,Sin,clock, sec_or_min, cell_rec, child, config, createRec, threshold, same_thresh) 	
	% config: config.adder decide adder or sizer, config.reset_S decide if S reset to 0 after cell division, config.singOper decide nature of Operon - delayed(1), synced(2), independent(0)
	% createRec flag: 0 - Don't create and pass the cell simulation records, only pass the cell birth/death times, div_durs.
	
	seed_store = rng;	% Only stores seed values, to check if any simulation iteration have same state or not
	% config
    adder = config.adder;
    reset_S = config.reset_S;
    singOper = config.singOper;
	
	% Load parameters
	%sec_or_min = 0; % sec
	%sec_or_min = 1; % min
	par_load_growthFN_v2;
	ncells = 1;
    
	% ES_ratio from S/(S+KM)
	S = cell(p,n);
	ES_ratio = zeros(p,length(0:0.05:1)-1);

	% Setup the initial values for the ODE
	x0 = cell_rec.x0;
	size_bir = x0(each_cell-p+(1:p));

	% Initialization of enz_profile variable
	cell_dat{1} = struct('burst_dat',{cell(p,n)},	'mrna_dat',{cell(p,n)},	'prot_dat', {cell(p,n)}, 'mrna_end',cell_rec.mrna_end,	'prot_end',cell_rec.prot_end, 'clock',clock, 'parent',cell_rec.parent, 'birth_time',cell_rec.div_t);
	
	% ODE settings
	opti1 = odeset('Events',@(t,x) cell_div_event(t,x,ncells,p,n,metab_threshold, each_cell),'NormControl','on','RelTol',1e-3);
	opti2 = odeset(opti1,'NonNegative',1:each_cell,'JPattern',JPattern_MM(ncells,n,p,0));	% No CrossFeeding
	
	% Loop to grow cells starting from 1 of each type (only 1 if 1 type), and grown till max_cell
	% clock = tspan(2);
	count = 0;
	divs = 0;
	div_t = clock;
	extend = 0;
	sol = struct('xe',[]);

	while isempty(sol.xe)
		% Create or Extend the profile
		cell_dat{1} = enz_prof_fn3(n,p,T,gene_seqlen,speed_transcription,t_on,t_off, prot_seqlen, speed_translation, clock, extend, cell_dat{1}, singOper);
		
		% Generate and Simulate the ODE function for metabolite production, and cross-feeding
		tspan = [clock clock+T];
		sol = ode15s(@(t,x) ODE_gen_MMfn_analysis_noFeed(t,x, 1, p, n, Sin, par, metab_threshold, cell_dat(1)), tspan, x0, opti2);
		count = count + 1;
		sol2{count} = sol;
		
		% Count the ES/ETot ratio
		eq_t = sol.x(1):(20/60):sol.x(end);		% every 20 seconds
		cell_id = 1;								% In case of single_cell_lineage = 1
		%cell_id = 1:each_cell:each_cell*ncells;	% In case of simulation of multiple_cells, need this
		for i = 1:ncells
			S = cell(p,n);
			for j = 1:p
				for k = 1:n
					S{j,k} = horzcat(S{j,k}, interp1(sol.x,sol.y(cell_id(i)+(k-1)*p+j-1,:), eq_t, 'linear', 'extrap'));
					ES_ratio(j,:) = ES_ratio(j,:) + histcounts(S{j,k}./(S{j,k} + KM),[0:0.05:1]);
				end
			end
		end
			
		clock = clock + T;
		x0 = sol.y(:,end);
		extend = 1;
	end
	
	div_t = sol.xe;
	% Post meeting division threshold
	% Remaining mrna and protein in the cell at division
	for i = 1:p
		for j = 1:n
			mrna_end(i,j) = interpl_lastval(cell_dat{1}.mrna_dat{i,j}.tot_rna_profile.t, cell_dat{1}.mrna_dat{i,j}.tot_rna_profile.v, div_t);
			prot_end(i,j) = interpl_lastval(cell_dat{1}.prot_dat{i,j}.tot_prot_profile.t, cell_dat{1}.prot_dat{i,j}.tot_prot_profile.v, div_t);
		end
	end
	
	% Create cell_rec for later analysis of data (createRec == 1 => Detailed record of profile. 0 => basic record)
	cell_rec = compile_cell_rec(n,p,cell_dat{1}.birth_time, div_t, cell_dat{1}, sol2, mrna_end, prot_end, ES_ratio, cell_rec.parent, cell_rec.mrna_end,cell_rec.prot_end, createRec);
	bac = cell_rec;
	
	% Initial parameters for the cell's 2 daughter cells
	bac.progeny.mrna = zeros(p,n,2);
	bac.progeny.prot = zeros(p,n,2);
	bac.progeny.x0 = zeros(each_cell,2);
	
	% Decide Cell Partition Ratio. Considering the cell is well-mixed all the contents are segregated by the given ratio.
	part_ratio = binornd(1000,0.5)/1000;
	
	% mrna distribution
	bac.progeny.mrna(:,:,1) = round(part_ratio*cell_rec.mrna_end);
	bac.progeny.mrna(:,:,2) = cell_rec.mrna_end - bac.progeny.mrna(:,:,1);
	
	% prot distribution
	bac.progeny.prot(:,:,1) = round(part_ratio*cell_rec.prot_end);
	bac.progeny.prot(:,:,2) = cell_rec.prot_end - bac.progeny.prot(:,:,1);
	
	% metabolite distribution
	if adder == 1
		% Pure Adder, reset final metabolites to ZERO
		x0(each_cell-p+(1:p),:) = 0;
	elseif adder == 0
		% Metabolic Sizer, metab_threshold amount gets used up, rest is passed on to the daughter cells
		x0(each_cell-p+(1:p)) = x0(each_cell-p+(1:p)) - metab_threshold(:);
	% New Sizer and Adder definitions. Do not require consumption of the metabolites
	%elseif adder == 2 |	adder == 3 
		%x0(each_cell-p+(1:p)) = x0(each_cell-p+(1:p));
	end
	if reset_S == 1		
		% Do not carry over the Substrate to daughter cells
		x0(1) = 0;		
	end
	bac.progeny.x0(:,1) = part_ratio*x0;
	bac.progeny.x0(:,2) = x0 - bac.progeny.x0(:,1);
	
	bac.seed_store = seed_store;
	bac.size_div = x0(each_cell-p+(1:p));
	bac.size_bir = size_bir;
end