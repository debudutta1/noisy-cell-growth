% Parameters for simulation

gene_seqlen = 1029; % nucleotides
prot_seqlen = 343; % amino acids

k3 = 10;			% per sec
KM = 60220;		% Conc units time independent

if sec_or_min == 1	% min
	speed_transcription = 55*60; 	% nucleotides/min
	speed_translation = 20*60; 		% amino acids/second
	
	% Exponential distribution D1 to derive the duration of transcription ON time
	% On Average transcription ON time in minutes
	genExp.t_ON = genExp.t_ON*1; 
	genExp.t_ON_ov = genExp.t_ON_ov*1; 
	% Exponential distribution D2 to derive the duration of transcription OFF time
	% On Average transcription OFF time in minutes
	genExp.t_OFF = genExp.t_OFF*1; 
	genExp.t_OFF_ov = genExp.t_OFF_ov*1; 
	
	T = T*60;   % Size of the enzyme time profile, in mins
	
	% Parameters for Enzyme Kinetics
	Sin = Sin*60;	% per minute 
	Pin = 60*1000*ones(n,p); 
	k3 = k3*60;			% per min
	
	par.feed_rate = feed_rate;	% In Noisy Transport feed_rate is a fraction
%	secrete_rate = secrete_rate*60;

elseif sec_or_min == 0	% sec
	speed_transcription = 55; 	% nucleotides/second	
	speed_translation = 20; 	% amino acids/second
	% Exponential distribution D1 to derive the duration of transcription ON time
	% On Average transcription ON time in minutes, converted to seconds
	genExp.t_ON = genExp.t_ON*60; 
	genExp.t_ON_ov = genExp.t_ON_ov*60; 
	% Exponential distribution D2 to derive the duration of transcription OFF time
	% On Average transcription OFF time in minutes, converted to seconds
	genExp.t_OFF = genExp.t_OFF*60; 
	genExp.t_OFF_ov = genExp.t_OFF_ov*60; 
	
	T = T*60*60;   % Size of the enzyme time profile, in seconds
	
	% Parameters for Enzyme Kinetics
	Sin = Sin;	% per sec 
	Pin = 1000*ones(n,p);
	k3 = k3;		% per sec
	
	par.feed_rate = feed_rate;	% In Noisy Transport feed_rate is a fraction
%	secrete_rate = secrete_rate;
end
	
% time to transcribe the gene
transcribe_time = gene_seqlen/speed_transcription;
% time to translate the transcript
translate_time = prot_seqlen/speed_translation;

par.h = 10; 
par.sec_KM = 1e5;
%par.f1 = 1/10000; 
par.f2 = 100;

% Secretion kinetics
par.secRatio = secRatio;
%par.k4 = repmat(secRatio,[1,2]);
par.KM = repmat(KM,p,n);


% Two cell types, complementary auxotrophs
par.k3 = repmat(k3,[2,p,n]);

% Decide the scheme of production of secreted metabolite: 1) Only Flux Doubling, 2) Only enzyme Overexpression, 3) Both 
if OVprod_mode == 1	% 1) Only Flux Doubling in Secreted pathway
	par.k3(2,1,:) = 2*par.k3(2,1,:);
	par.k3(1,2,:) = 2*par.k3(1,2,:);
elseif OVprod_mode == 2		% 2) Only enzyme Overexpression,
	% Two cell types, complementary auxotrophs
	par.k3 = repmat(k3,[2,p,n]);	
elseif OVprod_mode == 3		% 3) Both 
	par.k3(2,1,:) = 2*par.k3(2,1,:);
	par.k3(1,2,:) = 2*par.k3(1,2,:);
end
	
if exist('p0')	% Implies auxotrophy simulations, need p+1 metab_threshold
	p = p0;		% Temporary change to reduce code redundancy
end

metab_thres = threshold;
if same_thresh == 1
	if adder == 2	% New SIZER
		for i = 1:p
			metab_threshold(i) = metab_thres*2;		% New SIZER is twice that of earlier value since, now we do not consider that the produced value is used up. The starting values of metabolites represent a size proxy
		end
	elseif adder == 3	% New ADDER
		for i = 1:p
			metab_threshold(i) = cell_rec.x0(end-p+i) + metab_thres;
		end
	end
elseif same_thresh == 0
	if adder == 2	% New SIZER
		for i = 1:p
			metab_threshold(i) = threshold(i)*2;
		end
	elseif adder == 3	% New ADDER
		for i = 1:p
			metab_threshold(i) = cell_rec.x0(end-p+i) + threshold(i);
		end
	end
	for i = 1:p
		metab_threshold(i) = threshold(i);
	end
end

if exist('p0')	% Implies auxotrophy simulations, need p+1 metab_threshold
	p = p0 - 1;	% Revert to original value
	% ODE specifications size
	each_cell = (p*n+1+1);
else
	% ODE specifications size
	each_cell = (p*n+1);
end