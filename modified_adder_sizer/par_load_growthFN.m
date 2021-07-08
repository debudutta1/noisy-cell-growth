% Parameters for simulation
%sec_or_min = 1; % min
%sec_or_min = 0; % sec

gene_seqlen = 1305; % nucleotides
prot_seqlen = 434; % amino acids
% variables declared
%n = 3;   % Size of enzyme cascade
%p = 1;
%T = 0.25;   % In hrs 
%Sin = 11000*p;	% per sec 
%Sin = 19000*p;	% per sec 

if sec_or_min == 1	% min
	speed_transcription = 55*60; 	% nucleotides/min
	speed_translation = 20*60; 		% amino acids/min
	
	% Exponential distribution D1 to derive the duration of transcription ON time
	t_on = 6; % On Average transcription ON time in minutes
	% Exponential distribution D2 to derive the duration of transcription OFF time
	t_off = 37; % On Average transcription OFF time in minutes
	
	T = T*60;   % Size of the enzyme time profile, in mins
	
	% Parameters for Enzyme Kinetics
	Sin = Sin*60;	% per minute 
	%Pin = 60*1000*ones(n,p); 
	k3 = 12*60;			% per min
	
	load('auxPin.mat');
	par.auxPin = auxPin*60;

elseif sec_or_min == 0	% sec
	speed_transcription = 55; 	% nucleotides/second	
	speed_translation = 20; 	% amino acids/second
	% Exponential distribution D1 to derive the duration of transcription ON time
	t_on = 6*60; % On Average transcription ON time in minutes, converted to seconds
	% Exponential distribution D2 to derive the duration of transcription OFF time
	t_off = 37*60; % On Average transcription OFF time in minutes, converted to seconds
	
	T = T*60*60;   % Size of the enzyme time profile, in seconds
	
	% Parameters for Enzyme Kinetics
	Sin = Sin;	% per sec 
	%Pin = 1000*ones(n,p);
	k3 = 12;		% per sec
	
	load('auxPin.mat');
	par.auxPin = auxPin;
end
	
% time to transcribe the gene
transcribe_time = gene_seqlen/speed_transcription;
% time to translate the transcript
translate_time = prot_seqlen/speed_translation;

h = 10; 
f1 = 1/3; 
f2 = 100;

KM = 34347.826;		% Conc units time independent

par.k3 = repmat(k3,p,n);
par.k4 = repmat(Sin/10,p,n);
par.KM = repmat(KM,p,n);

% Auxotrophy
%par.auxPin = Sin/(p+1);
%par.auxPin = 20000;
%global auxPin;

	
% For gentime = 30 mins, metab_threshold =  1.9745e+07, Sin = 10969 ~ 11000
%Sin = 11000; % Influx rate of Substrate into the cell 
metab_thres = 1.9745e+07;
%metab_thres = 5.64e7;
%metab_thres = 1e+04;

if exist('p0')	% Implies auxotrophy simulations, need p+1 metab_threshold
	for i = 1:p0
		metab_threshold(i) = metab_thres;
	end
	% ODE specifications size
	each_cell = (p*n+1+1);
else
	for i = 1:p
		metab_threshold(i) = metab_thres;
	end
	% ODE specifications size
	each_cell = (p*n+1);
end

%	cf_size = 3*p*n;	% Old all CF size