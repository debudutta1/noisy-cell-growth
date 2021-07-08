% Lineage mapped cell division simulation
function y = parallel_growth_sim_ALT_auxFeedSec_noisyTransport(gens, Sin, n,p,T, config, sec_or_min, createRec, threshold, same_thresh, genExp, n_ext, feed_rate, cell_type, secRatio, OVprod_mode)

	% Load Single Lineage Simulation data to start the Exponential Simulation from
	%load('sinLin1000_start.mat');
	% OR SIMULATE 100 protocells!
	proto_n = 100;
	clock = 0;

	p1 = p - 1;
	each_cell = ((n*p1+1)+1) + 1; % For Secretion quantification
	Sin = Sin*p;	% Sin PASSED IS FOR PER CASCADE!
	
	% config = struct('adder',1,'reset_S',0,singOper,1);
	% Dummy ZERO Initialization for simulation start
	zero_cell_dat = struct('mrna_end',zeros(p1,n),'prot_end',zeros(p1,n),'parent',0,'birth_time',clock, 'x0', zeros(1, each_cell), 'index', 0, 'div_t', clock, 'type', cell_type,'T_mrna_end',0,'T_prot_end',0);
	child = 1;
	proto_cell{1} = cell_growth_fn_ALT_auxFeedSec_noisyTransport(n,p,T,Sin, clock, sec_or_min, zero_cell_dat, child, config, createRec, threshold, same_thresh, genExp, n_ext, feed_rate, cell_type, secRatio, OVprod_mode) ;
	y.div_durs(1) = proto_cell{1}.div_t-proto_cell{1}.birth_t;
	y.size_div(1,:) = proto_cell{1}.size_div;
	y.size_bir(1,:) = proto_cell{1}.size_bir;
	for i = 2:proto_n
		clock = proto_cell{i-1}.div_t;
		
		% Create input for cell of next generation
		next_dat = struct('mrna_end',proto_cell{i-1}.progeny.mrna(:,:,1),	'prot_end',proto_cell{i-1}.progeny.prot(:,:,1), 'parent',i-1, 'x0', proto_cell{i-1}.progeny.x0(:,1), 'index', i-1, 'div_t', proto_cell{i-1}.div_t, 'T_mrna_end',proto_cell{i-1}.progeny.T_mrna(1),	'T_prot_end',proto_cell{i-1}.progeny.T_prot(1));
		
		proto_cell{i} = cell_growth_fn_ALT_auxFeedSec_noisyTransport(n,p,T,Sin, clock, sec_or_min, next_dat, child, config, createRec, threshold, same_thresh, genExp, n_ext, feed_rate, cell_type, secRatio, OVprod_mode) ;
		y.div_durs(i) = proto_cell{i}.div_t-proto_cell{i}.birth_t;
		y.size_div(i,:) = proto_cell{i}.size_div;
		y.size_bir(i,:) = proto_cell{i}.size_bir;
		
		%display(y.div_durs(i));
		%disp(i);
	end
	
	% Adjust for optimal T!
	if mean(y.div_durs)/3600 > 0.8*T  
		T = 1.25*mean(y.div_durs)/3600;
	end
	
	% Exponential Growth Simulation results Record
	clock = proto_cell{proto_n}.div_t;
	% Create input for cell of next generation
	next_dat = struct('mrna_end',proto_cell{proto_n}.progeny.mrna(:,:,1),	'prot_end',proto_cell{proto_n}.progeny.prot(:,:,1), 'parent',i-1, 'x0', proto_cell{proto_n}.progeny.x0(:,1), 'index', i-1, 'div_t', proto_cell{proto_n}.div_t, 'T_mrna_end',proto_cell{proto_n}.progeny.T_mrna(1),	'T_prot_end',proto_cell{proto_n}.progeny.T_prot(1));
		
	cell_rec{1} = cell_growth_fn_ALT_auxFeedSec_noisyTransport(n,p,T,Sin, clock, sec_or_min, next_dat, child, config, createRec, threshold, same_thresh, genExp, n_ext, feed_rate, cell_type, secRatio, OVprod_mode) ;
	
	y.div_durs_exp(1) = cell_rec{1}.div_t - cell_rec{1}.birth_t;
	y.size_bir_exp(1,:) = cell_rec{1}.size_bir;
	y.size_div_exp(1,:) = cell_rec{1}.size_div;
	
	y.secretion_prof{1} = cell_rec{1}.secretion_prof;	% Take first data entry as cell array
	y.sec_end(1) = cell_rec{1}.sec_end;
	y.all_metab(:,1) = cell_rec{1}.all_metab;
	y.mrna_produced(:,:,1) = cell_rec{1}.mrna_produced;
	y.prot_produced(:,:,1) = cell_rec{1}.prot_produced;
	y.mrna_beg(:,:,1) = cell_rec{1}.mrna_beg;
	y.mrna_end(:,:,1) = cell_rec{1}.mrna_end;
	y.prot_beg(:,:,1) = cell_rec{1}.prot_beg;
	y.prot_end(:,:,1) = cell_rec{1}.prot_end;
	
	for i = 1:(gens-1)
		parfor j = 1:(2^i)
			index = 2^i - 1 + j;
			% Given the index of the current cell, calculate the index of the parent cell. Also determine which is the present index the 1st or 2nd progeny of the parent cell
			if rem(index,2) == 0	% even
				parent = index/2;	
				progeny = 1;
			else
				parent = (index - 1)/2;
				progeny = 2;
			end
			% Create input for cell of next generation
			next_dat = struct('mrna_end',cell_rec{parent}.progeny.mrna(:,:,progeny),	'prot_end',cell_rec{parent}.progeny.prot(:,:,progeny), 'parent',parent, 'x0', cell_rec{parent}.progeny.x0(:,progeny), 'index', index, 'div_t', cell_rec{parent}.div_t, 'T_mrna_end',cell_rec{parent}.progeny.T_mrna(progeny),	'T_prot_end',cell_rec{parent}.progeny.T_prot(progeny));
			
			clock = cell_rec{parent}.div_t;
			new_cell_rec{j} = cell_growth_fn_ALT_auxFeedSec_noisyTransport(n,p,T,Sin, clock, sec_or_min, next_dat, child, config, createRec, threshold, same_thresh, genExp, n_ext, feed_rate, cell_type, secRatio, OVprod_mode) ;
			
			div_durs_exp(j) = new_cell_rec{j}.div_t - new_cell_rec{j}.birth_t;
			size_bir_exp(j,:) = new_cell_rec{j}.size_bir;
			size_div_exp(j,:) = new_cell_rec{j}.size_div;
			
			secretion_prof{j} = new_cell_rec{j}.secretion_prof;
			sec_end(j) = new_cell_rec{j}.sec_end;
			all_metab(:,j) = new_cell_rec{j}.all_metab;
			mrna_produced(:,:,j) = new_cell_rec{j}.mrna_produced;
			prot_produced(:,:,j) = new_cell_rec{j}.prot_produced;
			mrna_beg(:,:,j) = new_cell_rec{j}.mrna_beg;
			mrna_end(:,:,j) = new_cell_rec{j}.mrna_end;
			prot_beg(:,:,j) = new_cell_rec{j}.prot_beg;
			prot_end(:,:,j) = new_cell_rec{j}.prot_end;
					
			%display(index);
		end
		cell_rec(2^i:2^(i+1)-1) = new_cell_rec(1:2^i);
		y.div_durs_exp(2^i:2^(i+1)-1) = div_durs_exp(1:2^i);
		y.size_bir_exp(2^i:2^(i+1)-1,:) = size_bir_exp(1:2^i,:);
		y.size_div_exp(2^i:2^(i+1)-1,:) = size_div_exp(1:2^i,:);
		
		y.secretion_prof(2^i:2^(i+1)-1) = secretion_prof(1:2^i);	% Transfer cell array as the cell array
		y.sec_end(2^i:2^(i+1)-1) = sec_end(1:2^i);
		y.all_metab(:,2^i:2^(i+1)-1) = all_metab(:,1:2^i);
		y.mrna_produced(:,:,2^i:2^(i+1)-1) = mrna_produced(:,:,1:2^i);
		y.prot_produced(:,:,2^i:2^(i+1)-1) = prot_produced(:,:,1:2^i);
		y.mrna_beg(:,:,2^i:2^(i+1)-1) = mrna_beg(:,:,1:2^i);
		y.mrna_end(:,:,2^i:2^(i+1)-1) = mrna_end(:,:,1:2^i);
		y.prot_beg(:,:,2^i:2^(i+1)-1) = prot_beg(:,:,1:2^i);
		y.prot_end(:,:,2^i:2^(i+1)-1) = prot_end(:,:,1:2^i);
			
		clear new_cell_rec;
		display(2^(i+1));
	end
	y.proto_cell = proto_cell;
	y.cell_rec = cell_rec;
end