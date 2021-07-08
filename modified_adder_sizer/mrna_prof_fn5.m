% Function to create a stochastic mRNA profile, by taking input of n - cascasde steps, T - time duration of profile, gene_seqlen -gene seq length, speed_transcription
% Function assumes that all the enzymes are on the same operon and thus creates a single stochastic profile and assumes that the subsequent enzymes expressed the same as the first but with a time delay equal to the size of the first gene transcribed.
function y = mrna_prof_fn5(T,mrna_life,gene_seqlen,speed_transcription, clock, mrna_dat, mrna_end, t1)%, prev_t, t_offset) 	
	% t1 - input burst profile, prev_t - previous burst profile, clock- start position of profile, 
	% mrna_dat - previous mrna profile, mrna_end - value of mrna at end of last profile which is being extended
	
	% Exponential distribution D3 to derive the duration for the decay of the transcript, which in turn determines the number of proteins made.
	%mrna_life = 25*60; % Average life time of a transcript in minutes, converted to seconds
	% time to transcribe the gene
	transcribe_time = gene_seqlen/speed_transcription;
	
	% mRNA Production profile
	%% Calculate given the time length T, a list of mRNA creation events and a list of decay events
	mrna_birth = []; mrna_death = [];
	
	% Since mrna Life-times are exponential, it does not matter if we randomly compute them at any point in time! In case of a new cell or extending the simulation of a cell!
	% Initialize with inherited RNA
	if mrna_end > 0
		if isempty(mrna_dat)	% Implies that the function is called in extend = 0 mode
			mrna_birth(1:mrna_end) = clock;
			mrna_death(1:mrna_end) = clock + exprnd(mrna_life, [1, mrna_end]);
		else				% Implies that the function is called in extend = 1 mode
			in = (mrna_dat.mrna_birth <= clock & mrna_dat.mrna_death > clock);
			if sum(in) == mrna_end
				mrna_birth(1:mrna_end) = mrna_dat.mrna_birth(in);
				mrna_death(1:mrna_end) = mrna_dat.mrna_death(in);
			else 
				error('Error: mrna_end does not match!')
			end
		end
	end
		
	%% In case, the profile calculated cut across a burst then the mrna that were produced in that after previous clock+T have to be accounted for!
	if ~isempty(mrna_dat)
		in = (mrna_dat.mrna_birth > clock);
		mrna_birth(end+(1:sum(in))) = mrna_dat.mrna_birth(in);
		mrna_death(end+(1:sum(in))) = mrna_dat.mrna_death(in);
	end
	
	%% SINCE we do incorporate the products of previous burst in extended profile, we must offset the burst profile to start from the end of last burst!
	% HOWEVER, this step is already taken care of in the calling function of the burst
	%% ALTERNATELY, CALL Burst gen function in this function and create accordingly, in which case t1 will be the previous burst profile!
	%t1 = burst_prof_fn6(prev_t.burst_end + t_offset, n, T, transcribe_time, t_on,t_off);
	
	% Iterate through each burst
	%for bur = 1:length(t1.start)
	for bur = find(t1.start >= clock & t1.start < clock+T)		% ITERATE THROUGH ONLY THE NEW BURST PROFILES that came after clock, since the truncated bursts are taken care of!
		event = t1.start(bur);	% Start from the first burst event
		%mrna = floor((t1.stop(bur) - t1.start(bur))/transcribe_time); % Total Number of mRNA molecules produced in the transcriptional burst
		m_births = (event+transcribe_time):transcribe_time:t1.stop(bur);	% Only produce proteins within the burst
		mrna = length(m_births); 	% Total Number of mRNA molecules produced in the transcriptional burst
		
		end_pos = length(mrna_death);
		if mrna > 0		% To avoid any error in case no mRNA are produced in this burst
			% birth and death times of mRNA
			%mrna_birth(end_pos+(1:mrna)) = (event+transcribe_time):transcribe_time:t1.stop(bur);	% Only produce proteins within the burst
			mrna_birth(end_pos+(1:mrna)) = m_births;
			mrna_death(end_pos+(1:mrna)) = mrna_birth(end_pos+(1:mrna)) + exprnd(mrna_life, [1 mrna]);
		end
	end
			
	% Creating the tot_rna_profile
	if ~isempty(mrna_birth)	% length of birth and death are identical any one works
		tot_rna_profile.t(1:length(mrna_birth)) = mrna_birth;
		tot_rna_profile.c(1:length(mrna_birth)) = +1;
		tot_rna_profile.t(end+(1:length(mrna_birth))) = mrna_death;
		tot_rna_profile.c(end+(1:length(mrna_birth))) = -1;
		
		% Create sorted profile
		[A, B] = sort(tot_rna_profile.t,'ComparisonMethod','real');
		% Limit the profile generated to clock+T or last mrna trancribed
		in = (A >= clock & A <= clock+T);
		A = A(in); 	
		B = B(in);
		tot_rna_profile.t = A;
		tot_rna_profile.c = tot_rna_profile.c(B);
		if sum(A == clock) == 0		% interp1 function also requires that the values are all unique, thus the data series must not have two repeating values!
			% Add the data point at clock in order to allow for interpolation function with 'last' to work in ODE
			tot_rna_profile.t = cat(2, clock, tot_rna_profile.t);
			tot_rna_profile.c = cat(2, mrna_end, tot_rna_profile.c);
		end	
	else%if isempty(mrna_birth)
		tot_rna_profile.t = clock;
		tot_rna_profile.c = mrna_end;
	end
	
	tot_rna_profile.v = cumsum(tot_rna_profile.c);
			
	% Additional closing entry to avoid interp1 error only 1 entry!
	if length(tot_rna_profile.t) == 1
		tot_rna_profile.t(end+1) = clock+T;
		tot_rna_profile.v(end+1) = tot_rna_profile.v(1);
	end	
	
	% Update cum_rna_profile
	if ~isempty(mrna_dat)
	% So that the end value of previous cum_rna_profile and the first value of the new profile are not repeated, since they are redundant, and non-uniqueness causes issues with 
		cum_rna_profile.t = cat(2, mrna_dat.cum_rna_profile.t, tot_rna_profile.t(2:end));	
		cum_rna_profile.v = cat(2, mrna_dat.cum_rna_profile.v, tot_rna_profile.v(2:end));	
	else
		cum_rna_profile.t = tot_rna_profile.t;	
		cum_rna_profile.v = tot_rna_profile.v;	
	end
	
	% Clear c filed from structure
	tot_rna_profile = rmfield(tot_rna_profile,'c');
	
	if tot_rna_profile.v(end) ~= sum(mrna_birth <= clock+T & mrna_death > clock+T)
		error('mrna_end doesnot match tot_rna_profile and mrna_birth/death')
	end
	y.mrna_end = tot_rna_profile.v(end);
	y.mrna_birth = mrna_birth;
	y.mrna_death = mrna_death;
	y.tot_rna_profile = tot_rna_profile;
	y.cum_rna_profile = cum_rna_profile;
end