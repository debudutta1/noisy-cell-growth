% Function to create BOTH the stochastic mRNA profile, and the corresponding Enzyme protein profile
% INPUT: n - cascasde steps, T - minimum time duration of profile, gene_seqlen, speed_transcription, t_ON, t_OFF, prot_seqlen, speed_translation, clock - starting time of the profile to generate, 

% Function creates a single stochastic profile based on input options singOper. fv4 is for CF simulations, where cell_dat.type variable exists additionally
function var = enz_prof_fn6_OV(n,p,T,gene_seqlen,speed_transcription, genExp, prot_seqlen, speed_translation, clock, extend, prev_cell_dat, singOper, OVprod_mode)
	% extend - if extending previous profile
	
	% Exponential distribution D3 to derive the duration for the decay of the transcript, which in turn determines the number of proteins made.
%	mrna_life = 25*60; 		% Average life time of a transcript in minutes, converted to seconds
	mrna_life = 2.5*60; 		% Average life time of a transcript in minutes, converted to seconds
	transcribe_time = gene_seqlen/speed_transcription; 	% time to transcribe the gene
	% % Exponential distribution D4 to derive the duration for the decay of the enzyme protein molecules
	prot_life = 20*60*60; % Average life time of a protein in minutes, converted to seconds
	translate_time = prot_seqlen/speed_translation; 	% time to translate the transcript
	
%	singOper = 2;	% 1 implies single operon, 2 implies synced, 0 implies independent expression for each of the n steps in cascade
	
	% Initialize cell_dat
	cell_dat =  prev_cell_dat;
	
	% For each of the p parallel cascades
	for i = 1:p
		% For each of the n serial cascades
		for j = 1:n
			
			% Decide the mRNA burst profile beginning
			if extend == 0
				burst_begin = clock;
				mrna_end = prev_cell_dat.mrna_end(i,j);
				prot_end = prev_cell_dat.prot_end(i,j);
			elseif extend == 1	
				if prev_cell_dat.burst_dat{i,j}.burst_end > clock
                    burst_begin = prev_cell_dat.burst_dat{i,j}.burst_end;
                else
                    burst_begin = clock;
                end
				mrna_end = prev_cell_dat.mrna_dat{i,j}.mrna_end;
				prot_end = prev_cell_dat.prot_dat{i,j}.prot_end;
			end
			
			% Decide OverProduction fom OVprod_mode
			% if OVprod_mode == 1
			% 
				t_on = genExp.t_ON;
				t_off = genExp.t_OFF;
			if OVprod_mode == 2 || OVprod_mode == 3
				% Index of secreted and auxotrophy
				sec_p = prev_cell_dat.type(3);
				if i == sec_p
					t_on = genExp.t_ON_ov;
					t_off = genExp.t_OFF_ov;
				end
			end
			
			% Decide the mRNA burst profile,
			if singOper == 0		% Independent gene expression 
				t = burst_prof_fn6(burst_begin, n, clock+T-burst_begin, transcribe_time, t_on,t_off);
				cell_dat.burst_dat{i,j} = t;
			elseif singOper == 1	% Single Operon with burst delays
				if j == 1
					t = burst_prof_fn6(burst_begin, n, clock+T-burst_begin, transcribe_time, t_on,t_off);
				end
				% OFFSET for each of the cascade steps
				t1.start = t.start + transcribe_time*(j-1);
				t1.stop = t.stop + transcribe_time*(j-1);
				t1.burst_end = t.burst_end + transcribe_time*(j-1);
				cell_dat.burst_dat{i,j} = t1;
			elseif singOper == 2	% All synced bursts
				if j == 1
					t = burst_prof_fn6(burst_begin, n, clock+T-burst_begin, transcribe_time, t_on,t_off);
                end
                cell_dat.burst_dat{i,j} = t;
			end
			
			%=======
			
			%% NEED TO CAREFULLY CREATE PROPER cell_dat for input in case of extend = 0 or 1, in CALLING FN
			cell_dat.mrna_dat{i,j} = mrna_prof_fn5(T, mrna_life, gene_seqlen, speed_transcription, clock, prev_cell_dat.mrna_dat{i,j}, mrna_end, cell_dat.burst_dat{i,j});
			
			cell_dat.prot_dat{i,j} = prot_prof_fn4(T, prot_life, prot_seqlen, speed_translation, clock, prev_cell_dat.prot_dat{i,j}, prot_end, cell_dat.mrna_dat{i,j});
		end
	end
	% Update clock
	cell_dat.clock = clock+T;
	
	% Variables to be kept constant and passed on till next cell division
	cell_dat.parent = prev_cell_dat.parent;
	cell_dat.type = prev_cell_dat.type;
	if extend == 0
		cell_dat.birth_time = clock;
%	else
%		cell_dat.birth_time = prev_cell_dat.birth_time;
	end
	
    % Set the start time of the cell birth
%	if extend == 0
%		y.birth = clock;
%	elseif extend == 1
%		y.birth = cell_dat.birth;
%	end
	
	%%%%%==========================
	var = cell_dat;
	
end