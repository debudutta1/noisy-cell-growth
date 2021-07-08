% MONOD'S OBSERVATION - VARY n
% RUN SIMULATION



p = [1 2 3];
Sin = 1e3*[5 6 6.5 7 9 13];
n = [1 2 3 5 10];

% define seed_range before running script
if exist('seed_range') == 0
	seed_range = 0:3;	
end

gens = 13; % 2^12 cells = 8192 ; 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
%n = 3;

T = 1.5;	% initial T in hrs

createRec = 0;
%createRec = 2;	% Detailed data not output, but computed to extract important parameters

same_thresh = 1;
if same_thresh == 1
	threshold = 1e+07;
else
	for i = 1:p
		threshold(i) = 1e+07;
	end
end

base_dir = 'D:\Debu Simulations\Sep 2020\var_n\';


% SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
fold_name = 'Sizer';


% for n == 1
for k = seed_range+1
	for l = 1:length(n)
		for j = 1%1:length(p)
			for i = length(Sin):-1:1
				
				load(strcat(base_dir,fold_name,'\var_n',num2str(n(l)),'_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k-1),'.mat'),'div_durs_exp','sim_vars','config','size_bir_exp','size_div_exp');
				
				div_durs_compiled(:,k,j,l,i) = div_durs_exp;
				size_bir_compiled(:,k,j,l,i) = size_bir_exp;
				size_div_compiled(:,k,j,l,i) = size_div_exp;
				
			end
		end
	end
end
% Obtain Growth rate of cells from dataset
parpool(length(seed_range))
growth_rate = nan(length(seed_range),1,length(n),length(Sin));
reps = 100;
for l = 1:length(n)
%	for j = 1:length(p)
	for j = 1
		for i = 1:length(Sin)
			parfor k = seed_range+1
				tic; growth_rate(k,j,l,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,l,i)); toc;
			end
		end
	end
end


% For n = 3
for k = seed_range+1
	for l = 1:length(n)
		for j = 3%1:length(p)
			for i = length(Sin):-1:1
				
				load(strcat(base_dir,fold_name,'\var_n',num2str(n(l)),'_p',num2str(p(j)),'_Sin',num2str(Sin(i)),'_rng',num2str(k-1),'.mat'),'div_durs_exp','sim_vars','config','size_bir_exp','size_div_exp');
				
				div_durs_compiled(:,k,1,l,i) = div_durs_exp;
				size_bir_compiled(:,:,k,1,l,i) = size_bir_exp;
				size_div_compiled(:,:,k,1,l,i) = size_div_exp;
				
			end
		end
	end
end
% Obtain Growth rate of cells from dataset
%parpool(length(seed_range))
growth_rate = nan(length(seed_range),1,length(n),length(Sin));
reps = 100;
for l = 1:length(n)
%	for j = 1:length(p)
	for j = 1
		for i = 1:length(Sin)
			parfor k = seed_range+1
				tic; growth_rate(k,j,l,i) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,l,i)); toc;
			end
		end
	end
end

%save(strcat(base_dir,'var_n_dat.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate');

%%===========================
lin_colors = lines(14);

f = figure; hold on
%j = 1;	% => p = 1
j = 1;	% => p = 1
for l = 1:length(n)
	plot(Sin, reshape(mean(growth_rate(:,j,l,:),1),length(Sin),[]),'-o','Displayname',strjoin({'n =',num2str(n(l))},' '),'linewidth',2,'color',lin_colors(l,:));
end
legend('location','southeast','NumColumns',1);
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux (per sec)')

f.Position(3:4) = f.Position(3:4)*0.75;
f.Children(2).XLabel.FontSize = 11;
f.Children(2).YLabel.FontSize = 11;

%name = 'vary_n_growR_sin_p1';
name = 'vary_n_growR_sin_p3';
set(f, 'Color', 'none');
export_fig(name, '-png', '-pdf','-r300', '-transparent', '-painters')


%%========================

f = figure; hold on
j = 1;	% => p = 1
for i = 1:length(Sin)
	plot(n, reshape(mean(growth_rate(:,j,:,i),1),length(n),[]),'-o','Displayname',strjoin({'Sin =',num2str(Sin(i))},' '),'linewidth',2,'color',lin_colors(i,:));
end
legend('location','northeast','NumColumns',1);
ylabel('Mean growth rate (per hour)')
xlabel('Number of enzyme steps in pathway')

% vary p
f = figure; hold on
j = 1;	% => p = 1
for j = 1:length(p)
	plot(Sin, reshape(mean(growth_rate(:,j,l,:),1),length(Sin),[]),'-o','Displayname',strjoin({'p =',num2str(p(j))},' '),'linewidth',2,'color',lin_colors(j,:));
end
legend('location','southeast','NumColumns',1);
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux (per sec)')