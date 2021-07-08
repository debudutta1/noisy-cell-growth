%% load data Fig 6a
sec_or_min = 0; % sec
n = 3;
p = 3; 
n_ext_range = 1;
feed_rate = [0 312 625 1250:1250:5000 7500 12500];
n_ext = 1;
k3_var_range = [0 1 2.5 5 10 25 50];
Sin = 1e3*[1.5 5 6 6.5 7 9 13];
seed_range = 0:1;

base_dir = 'D:\Debu Simulations\Sep 2020\var_extFeed_Rev1\';
load(strcat(base_dir,'var_extFeed_rev1_var_k3_withNoFeedback.mat'));

%%ANOVA	Fig 6a 
n_ext = 1;
%Sin_indices = 2:end;
%  div_durs_compiled(:,:,j,n_ext,i,k3_var)		=> j-> feed_rate, i -> Sin, 
j = 9; % feed = 12500
raw_dat = squeeze(div_durs_compiled(:,:,j,n_ext,2:end,:));

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

name = 'fig6a_anova';
writecell(tbl,strcat(name,'.csv'))

results = multcompare(stats,'Dimension',[1 2])