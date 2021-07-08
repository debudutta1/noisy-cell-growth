% Run ALTERNATIVE Auxotrophic Feed Simulations
% AIM: Estimate mean secretion rate that auxotrophic cells can sustain without taking a growth rate hit.
% APPROACH: 
% (1) FluxDouble: Set k_cat of pathway x2. Assume secR = 1/2 of all the produced final metabolite is secreted. (Tests without any change in noise profile (i.e. No effect on inter-division distribution, what is secretion rate possible), how much secretion capacity can be created by only flux doubling)
% (2) EnzDouble: Use 2 selected parameter sets: One where per cell protein expression ratio is 2:1, and another where the protein expression ratio is 2:1 calculated over the population. (Again by not taking a toll on the growth rate how much can be secreted?)
%% Vary secR = 0.3:0.1:0.7


sec_or_min = 0; % sec
n = 3;
p = 3;
n_ext = 0;

% define seed_range before running script
if exist('seed_range') == 0
	seed_range = 0:3;	
end

%Sin = 1e3*[1.5 4 5.5 6.5 7 7.5 10 20];
%Sin = 1e3*[1.5 3 4 5 6 6.5 7 9 13];
Sin = 1e3*[4 5 6 6.5 7 9];

% Noisy Feed rates are fractions to multiply with k3.
feed_rate = [0.5 0.8 0.9 1 1.1 1.2 2];

secRatio = 0.5;	% Only 0.5 to compare with the non-noisy CF case

OVprod_mode = 1; % Decide the scheme of production of secreted metabolite: 1) Only Flux Doubling, 2) Only enzyme Overexpression, 3) Both 

% Exponential distribution D1 to derive the duration of transcription ON time
t_ON = 4; % On Average transcription ON time in minutes
% Exponential distribution D2 to derive the duration of transcription OFF time
t_OFF = 2.4; % On Average transcription OFF time in minutes

genExp.t_ON = t_ON;
genExp.t_OFF = t_OFF;
% NO Overexpression
genExp.t_ON_ov = t_ON;
genExp.t_OFF_ov = t_OFF;

threshold = 1e7;
same_thresh = 1;

%addpath('D:\Google Drive\crossfeed_code\universal code_April2019\crossfeeding_ODE\cleaned final code\MM\Pure_Adder\New folder\modified_adder_sizer\');

% NEW ADDER with no delay - Synced
%config = struct('adder',3,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
% NEW SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent

gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
%createRec = 0;
createRec = 2;	% Detailed data not output, but computed to extract important parameters

% ALT version: Secretion from final metabolite 
cell_type1 = [1, 3, 2, 3];	%[1, 2, 2, 3];	%[type, sec_n, sec_p, aux_p]

% LAB WORKSTATION
addpath('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned')
% LAB PC
%addpath()

base_dir = 'D:\Debu Simulations\Sep 2020\';
mkdir(strcat(base_dir,'var_auxFeedSecrete_noisyTransport'));
base_dir = strcat(base_dir,'var_auxFeedSecrete_noisyTransport\');

%diary ALT_auxFeedSec_.log

for k = seed_range	%seed_range
	for i = 1:length(Sin)
		for j = length(feed_rate):-1:1
			for l = 1:length(secRatio)
				rng(k)
				
				display(strcat('seed_',num2str(k),'feed_',num2str(j),'sec_',num2str(l),'Sin_',num2str(i)))
				% If a simulation is already done and file exists, skip simulation for that condition
				if isfile(strcat(base_dir,'\ALT_aux1_Sin',num2str(Sin(i)),'_secR',num2str(secRatio(l)),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat'))
					display(strcat(base_dir,'\ALT_aux1_Sin',num2str(Sin(i)),'_secR',num2str(secRatio(l)),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat'))
					continue;
				end

				y = parallel_growth_sim_ALT_auxFeedSec_noisyTransport(gens, Sin(i), n, p, T, config, sec_or_min, createRec, threshold, same_thresh, genExp, n_ext, feed_rate(j), cell_type1, secRatio(l), OVprod_mode);

				div_durs = y.div_durs;
				div_durs_exp = y.div_durs_exp;
				proto_cell = y.proto_cell;
				cell_rec = y.cell_rec;
				size_bir = y.size_bir;
				size_div = y.size_div;
				size_bir_exp = y.size_bir_exp;
				size_div_exp = y.size_div_exp;
				
				sec_end = y.sec_end;
				secretion_prof = y.secretion_prof;
				all_metab = y.all_metab;
				mrna_produced = y.mrna_produced;
				prot_produced = y.prot_produced;
				mrna_beg = y.mrna_beg;
				mrna_end = y.mrna_end;
				prot_beg = y.prot_beg;
				prot_end = y.prot_end;
				
				sim_vars = struct('Sin',Sin(i),'p',p,'secRatio',secRatio(l),'feed_rate',feed_rate(j),'seed',k,'OVprod_mode', OVprod_mode, 'genExp', genExp);
					
				save(strcat(base_dir,'\ALT_aux1_Sin',num2str(Sin(i)),'_secR',num2str(secRatio(l)),'_feed',num2str(feed_rate(j)),'_rng',num2str(k),'.mat'),'div_durs_exp','sim_vars','config','size_bir_exp','size_div_exp','secretion_prof','all_metab','mrna_beg','mrna_end','mrna_produced','prot_beg','prot_end','prot_produced','sec_end');
			end
		end
	end
end

%diary off