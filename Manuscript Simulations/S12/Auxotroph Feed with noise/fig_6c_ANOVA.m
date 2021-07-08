%% load data Fig 6c
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

base_dir = 'D:\Debu Simulations\Sep 2020\';
base_dir = strcat(base_dir,'var_auxFeedSecrete_noisyTransport\');
load(strcat(base_dir,'rev2_auxFeedSec_noisy_dat.mat'),'div_durs_compiled');



%%ANOVA	Fig 6c 
plot_indices = [1:4 6];	% feed_ratio
sel_dat = 2:length(Sin);
%  div_durs_compiled(:,:,j,n_ext,i,k3_var)		=> j-> feed_rate, i -> Sin, 
j = 9; % feed = 12500
raw_dat = squeeze(div_durs_compiled(:,:,plot_indices,sel_dat));

dims = size(raw_dat);
%dims = dims([1 2 4 3]);
clear data
data = nan(prod(dims(1:3)),dims(4));

ct = 0;
for i = 1:dims(3)		
	for k = 1:dims(2)
		for l = 1:dims(1)
			ct = ct + 1;
			%data(ct,:) = squeeze(raw_dat(l,k,:,i));
			data(ct,:) = squeeze(raw_dat(l,k,i,:));
		end
	end
end


[p1,tbl,stats] = anova2(data,prod(dims(1:2)))

name = 'fig6c_anova';
writecell(tbl,strcat(name,'.csv'))

results = multcompare(stats,'Dimension',[1 2])