% Parameters for simulation
%sec_or_min = 1; % min
%sec_or_min = 0; % sec

%gene_seqlen = 1305; % nucleotides
%prot_seqlen = 434; % amino acids
gene_seqlen = 1029; % nucleotides
prot_seqlen = 343; % amino acids
% variables declared
%n = 3;   % Size of enzyme cascade
%p = 1;
%T = 0.25;   % In hrs 
%Sin = 11000*p;	% per sec 
%Sin = 19000*p;	% per sec 

% Exponential distribution D1 to derive the duration of transcription ON time
%t_on = 6; % On Average transcription ON time in minutes
% Exponential distribution D2 to derive the duration of transcription OFF time
%t_off = 6; % On Average transcription OFF time in minutes
%t_off = 3; % On Average transcription OFF time in minutes

% Enzyme parameters
%k3 = 12;			% per sec
%KM = 34347.826;		% Conc units time independent
k3 = 10;			% per sec
KM = 60220;		% Conc units time independent

% ODE specifications size
each_cell = (p*n+1);
cf_size = 3*p*n;

if sec_or_min == 1	% min
	speed_transcription = 55*60; 	% nucleotides/min
	speed_translation = 20*60; 		% amino acids/second
	
	% Exponential distribution D1 to derive the duration of transcription ON time
	t_on = t_on*1; % On Average transcription ON time in minutes
	% Exponential distribution D2 to derive the duration of transcription OFF time
	t_off = t_off*1; % On Average transcription OFF time in minutes
	
	T = T*60;   % Size of the enzyme time profile, in mins
	
	% Parameters for Enzyme Kinetics
	Sin = Sin*60;	% per minute 
	Pin = 60*1000*ones(n,p); 
	k3 = k3*60;			% per min
	
	if exist('feed_rate')
		par.feed_rate = feed_rate*60;
	end

elseif sec_or_min == 0	% sec
	speed_transcription = 55; 	% nucleotides/second	
	speed_translation = 20; 	% amino acids/second
	% Exponential distribution D1 to derive the duration of transcription ON time
	t_on = t_on*60; % On Average transcription ON time in minutes, converted to seconds
	% Exponential distribution D2 to derive the duration of transcription OFF time
	t_off = t_off*60; % On Average transcription OFF time in minutes, converted to seconds
	
	T = T*60*60;   % Size of the enzyme time profile, in seconds
	
	% Parameters for Enzyme Kinetics
	Sin = Sin;	% per sec 
	Pin = 1000*ones(n,p);
	k3 = k3;		% per sec
	
	if exist('feed_rate')
		par.feed_rate = feed_rate;
	end
end
	
% time to transcribe the gene
transcribe_time = gene_seqlen/speed_transcription;
% time to translate the transcript
translate_time = prot_seqlen/speed_translation;

h = 10; 
f1 = 1/3; 
f2 = 100;

par.k3 = repmat(k3,p,n);
par.k4 = repmat(Sin/10,p,n);
par.KM = repmat(KM,p,n);
	
% For gentime = 30 mins, metab_threshold =  1.9745e+07, Sin = 10969 ~ 11000
%Sin = 11000; % Influx rate of Substrate into the cell 
%metab_thres = 1.9745e+07;
%metab_thres = 5.64e7;
%metab_thres = 1e+04;

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

