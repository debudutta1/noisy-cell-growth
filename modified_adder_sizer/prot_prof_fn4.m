% Function to create a prot profile, by taking input of n - cascasde steps, T - time duration of profile, mrna profile
function y = prot_prof_fn4(T, prot_life, prot_seqlen,speed_translation,clock, prot_dat, prot_end, mrna_dat) % x1,x2 - prot_birth, prot_death, incom_prot - incomplete proteins formed before 	
	% clock- to locate position, x1 - prot_birth previous, x2 - prot_death previous
	
	% prot_seqlen = 434; % amino acids
	% speed_translation = 20; % amino acids/second
	translate_time = prot_seqlen/speed_translation; % time to translate the transcript
	% % Exponential distribution D4 to derive the duration for the decay of the enzyme protein molecules
	%prot_life = 20*60*60; % Average life time of a protein in minutes, converted to seconds
	% Protein Production
	%% Calculate given the time length T, a list of mRNA creation events and a list of decay events
	prot_birth = []; prot_death = [];
	
	% Since prot Life-times are exponential, it does not matter if we randomly compute them at any point in time! In case of a new cell or extending the simulation of a cell!
	% Initialize with inherited proteins
	if prot_end > 0
		if isempty(prot_dat)
			prot_birth(1:prot_end) = clock;
			prot_death(1:prot_end) = clock + exprnd(prot_life, [1, prot_end]);
		else
			in = (prot_dat.prot_birth <= clock & prot_dat.prot_death > clock);
			if sum(in) == prot_end
				prot_birth(1:prot_end) = prot_dat.prot_birth(in);
				prot_death(1:prot_end) = prot_dat.prot_death(in);
			else 
				error('Error: prot_end does not match!')
			end
		end
	end
	
	%% In case, the profile calculated cut across the lifetime of an mrna then the mrna that were produced in that after previous clock+T have to be accounted for!
	if ~isempty(prot_dat)
		in = (prot_dat.prot_birth >= clock);
		prot_birth(end+(1:sum(in))) = prot_dat.prot_birth(in);
		prot_death(end+(1:sum(in))) = prot_dat.prot_death(in);
	end

	% Iterate through each mrna
	%for bur = 1:length(mrna_dat.mrna_birth)		% ITERATE THROUGH ONLY THE NEW MRNA PROFILE
	for bur = find(mrna_dat.mrna_birth >= clock & mrna_dat.mrna_birth < clock+T)		% ITERATE THROUGH ONLY THE NEW MRNA PROFILE
		event = mrna_dat.mrna_birth(bur);	% Start from the first burst event
		%prot = floor((mrna_dat.mrna_death(bur) - mrna_dat.mrna_birth(bur))/translate_time); % Total Number of mRNA molecules produced in the transcriptional burst
		p_births = (event+translate_time):translate_time:mrna_dat.mrna_death(bur);
		prot = length(p_births); % Total Number of mRNA molecules produced in the transcriptional burst
		
		end_pos = length(prot_death);
		if prot > 0		% To avoid any error in case no mRNA are produced in this burst
			% birth and death times of mRNA
			%prot_birth(end_pos+(1:prot)) = (event+translate_time):translate_time:mrna_dat.mrna_death(bur);	% Only produce proteins within the burst
			prot_birth(end_pos+(1:prot)) = p_births;	% Only produce proteins within the burst
			prot_death(end_pos+(1:prot)) = prot_birth(end_pos+(1:prot)) + exprnd(prot_life, [1 prot]);
		end
	end
	
	% Creating the tot_prot_profile
	if ~isempty(prot_birth)	% length of birth and death are identical any one works		
		tot_prot_profile.t(1:length(prot_birth)) = prot_birth;
		tot_prot_profile.c(1:length(prot_birth)) = +1;
		tot_prot_profile.t(end+(1:length(prot_birth))) = prot_death;
		tot_prot_profile.c(end+(1:length(prot_birth))) = -1;
		
		% Create sorted profile
		[A, B] = sort(tot_prot_profile.t,'ComparisonMethod','real');
		% Limit the profile generated to clock+T or last prot trancribed
		in = (A >= clock & A <= clock+T);
		A = A(in); 	
		B = B(in);
		tot_prot_profile.t = A;
		tot_prot_profile.c = tot_prot_profile.c(B);
		if sum(A == clock) == 0		% interp1 function also requires that the values are all unique, thus the data series must not have two repeating values!
			% Add the data point at clock in order to allow for interpolation function with 'last' to work in ODE
			tot_prot_profile.t = cat(2, clock, tot_prot_profile.t);
			tot_prot_profile.c = cat(2, prot_end, tot_prot_profile.c);
		end		
	else%if isempty(mrna_birth)
		tot_prot_profile.t = clock;
		tot_prot_profile.c = prot_end;
	end
	
	tot_prot_profile.v = cumsum(tot_prot_profile.c);
	
	% Additional closing entry to avoid interp1 error only 1 entry!
	if length(tot_prot_profile.t) == 1
		tot_prot_profile.t(end+1) = clock+T;
		tot_prot_profile.v(end+1) = tot_prot_profile.v(1);
	end	
	
	% Update cum_prot_profile
	if ~isempty(prot_dat)
	% So that the end value of previous cum_prot_profile and the first value of the new profile are not repeated, since they are redundant, and non-uniqueness causes issues with 
		cum_prot_profile.t = cat(2, prot_dat.cum_prot_profile.t, tot_prot_profile.t(2:end));	
		cum_prot_profile.v = cat(2, prot_dat.cum_prot_profile.v, tot_prot_profile.v(2:end));	
	else
		cum_prot_profile.t = tot_prot_profile.t;	
		cum_prot_profile.v = tot_prot_profile.v;	
	end
	
	% Clear c filed from structure
	tot_prot_profile = rmfield(tot_prot_profile,'c');
	
	if tot_prot_profile.v(end) ~= sum(prot_birth <= clock+T & prot_death > clock+T)
		error('prot_end doesnot match tot_prot_profile and prot_birth/death')
	end
	y.prot_end = tot_prot_profile.v(end);
	y.prot_birth = prot_birth;
	y.prot_death = prot_death;
	y.tot_prot_profile = tot_prot_profile;
	y.cum_prot_profile = cum_prot_profile;
end