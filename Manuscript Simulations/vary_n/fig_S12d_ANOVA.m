%% load data Fig 6d
p = [1 2 3];
Sin = 1e3*[5 6 6.5 7 9 13];
n = [1 2 3 5 10];

% define seed_range before running script
if exist('seed_range') == 0
	seed_range = 0:3;	
end

base_dir = 'D:\Debu Simulations\Sep 2020\var_n\';
fold_name = 'Sizer';

% for n == 1
for k = seed_range+1
	for l = 1:length(n)
		for j = 1%1:length(p)
			for i = length(Sin):-1:1
				
				load(strcat(base_dir,fold_name,'\var_n',num2str(n(l)),'_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k-1),'.mat'),'div_durs_exp','sim_vars','config','size_bir_exp','size_div_exp');
				
				div_durs_compiled(:,k,l,i) = div_durs_exp;				
			end
		end
	end
end


%%ANOVA	Fig 6d 
raw_dat = squeeze(div_durs_compiled);

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

name = 'fig6d_anova';
writecell(tbl,strcat(name,'.csv'))

results = multcompare(stats,'Dimension',[1 2])