% Analysis code for ext Feed Simulations with varying k3 - Reviewer 1 suggestion

sec_or_min = 0; % sec
n = 3;
p = 3; 
n_ext_range = 1;
feed_rate = [0 312 625 1250:1250:5000 7500 12500];
n_ext = 1;
k3_var_range = [0 1 2.5 5 10 25 50];
Sin = 1e3*[1.5 5 6 6.5 7 9 13];
seed_range = 0:1;

% Exponential distribution D1 to derive the duration of transcription ON time
t_ON = 4; % On Average transcription ON time in minutes
% Exponential distribution D2 to derive the duration of transcription OFF time
t_OFF = 2.4; % On Average transcription OFF time in minutes

threshold = 1e+07;
same_thresh = 1;

addpath('C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned')

% NEW ADDER with no delay - Synced
%config = struct('adder',3,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent
% NEW SIZER with no delay - Synced
config = struct('adder',2,'reset_S',0,'singOper',2);	% singOper 2 synced, 1 delayed, 0 independent


gens = 13; % 2^14 cells = 16384
sec_or_min = 0; % 0 - sec
T = 1.5;
createRec = 0;

base_dir = 'D:\Debu Simulations\Sep 2020\var_extFeed_Rev1\';

for k3_var = 1:length(k3_var_range)
	for j = 1:length(feed_rate)
		for n_ext = n_ext_range
			% If no feed, n_ext doesn't matter. Skip redundant simulations
			if feed_rate(j) == 0 & n_ext > 1
				break;
			end
			for i = 1:length(Sin)
				for k = seed_range+1
	
					load(strcat(base_dir,'var_Sin',num2str(Sin(i)),'_nExt',num2str(n_ext),'_feed',num2str(feed_rate(j)),'_rng',num2str(k-1),'_k3feedbk_',num2str(k3_var_range(k3_var)),'.mat'),'div_durs_exp','size_bir_exp','size_div_exp');
					
					div_durs_compiled(:,k,j,n_ext,i,k3_var) = div_durs_exp;
					size_bir_compiled(:,:,k,j,n_ext,i,k3_var) = size_bir_exp;
					size_div_compiled(:,:,k,j,n_ext,i,k3_var) = size_div_exp;
					
				end
			end
		end
	end
end

% Obtain Growth rate of cells from dataset
parpool(3)
growth_rate = nan(length(seed_range),length(feed_rate),length(n_ext_range),length(Sin));
reps = 100;
for k3_var = 1:length(k3_var_range)
	for j = 1:length(feed_rate)
		for n_ext = n_ext_range
			if j == 1 && n_ext > 1
				break;
			end
			for i = 1:length(Sin)
				parfor k = seed_range+1
					tic; growth_rate(k,j,n_ext,i,k3_var) = exp_grow_rate(reps, gens, div_durs_compiled(:,k,j,n_ext,i,k3_var)); toc;
				end
			end
		end
	end
end

%save(strcat(base_dir,'var_extFeed_fig4.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate');
save(strcat(base_dir,'var_extFeed_rev1_var_k3.mat'),'div_durs_compiled','size_bir_compiled','size_div_compiled','growth_rate');

%=======================
base_dir = 'D:\Debu Simulations\Sep 2020\var_extFeed_Rev1\';
load(strcat(base_dir,'var_extFeed_rev1_var_k3_withNoFeedback.mat'));
%=======================


% Figure: Histogram of varying k3 for particular feed rate

lin_colors = lines(14);
% div_durs_compiled(:,k,j,n_ext,i,k3_var). i=> Sin, k => Seed, n_ext => No of metabs imported, j => feed_rate
n_ext = 1;

f = figure; 
hold on;
%i = 1; 	% Sin = 1.5e3
i = 2; % Sin = 5e3
j = 2;	% fixed feed rate = 312
for k3_var = 1:length(k3_var_range)
	c2(:,:) = div_durs_compiled(:,:,j,n_ext,i,k3_var)/60;
	
	[h1,h2] = histcounts(c2(:),'normalization','pdf');
	h1 = [0 h1];
	h2 = [h2(1) h2(1:end-1)+diff(h2)/2];	% Mean of the edges
	h(j) = fill(h2,h1,lin_colors(k3_var,:),'FaceAlpha',0.5,'EdgeColor','none','DisplayName',strjoin({'% feedback =',num2str(k3_var_range(k3_var))},' '));
end

%xlim([25 55])
xlabel('Cell generation times (in mins)')
ylabel('PDF')
%legend('location','northeast')



f = figure; 
hold on;
%i = 1; 	% Sin = 1.5e3
i = 2; % Sin = 5e3
j = 7;	% fixed feed rate = 312
clear('c2')
k = seed_range+1;
for k3_var = 1:length(k3_var_range)
	c2(:,:) = div_durs_compiled(:,k,j,n_ext,i,k3_var)/60;
	subplot(3,2,k3_var);
	histogram(c2(:),'normalization','pdf','DisplayName',strjoin({'% feedback =',num2str(k3_var_range(k3_var)),'Mean=',num2str(mean(c2(:)))},' '));
	legend('location','northeast')
end
xlabel('Cell generation times (in mins)')
ylabel('PDF')


%=======================
% Load data with no feedback
path_loc = 'C:\Users\sysbio_admin\Google Drive\EvoCoopSim_codes\Cell_Growth_Model_codeArchive\cleaned\Manuscript simulations\Monod\';
%ld_dat = load(strcat(base_dir,'var_extFeed_fig4_seed_0-6.mat'),'growth_rate','div_durs_compiled');%,'Sin','p','seed_range');
ld_dat = load('D:\Debu Simulations\Sep 2020\var_extFeed\var_extFeed_fig4_seed_0-6.mat','growth_rate','div_durs_compiled');%,'Sin','p','seed_range');
growth_rate_o = ld_dat.growth_rate;
clear('ld_dat');
load(strcat(path_loc,'divD_growR_vary_sin_p_sizer.mat'),'x_dat');%,'Sin','p','seed_range');
%=======================

% AuxFeed
auxfeed = load('D:\Debu Simulations\Sep 2020\var_auxFeed\auxFeed_fig5_dat2.mat','div_durs_compiled','growth_rate'); 
auxfeed.Sin = 1e3*[1.5 3 4 5 6 6.5 7 9 13];

%=======================
% Figure growth rate vs Sin

lin_colors = lines(14);

f = figure; hold on;

% Plot Auxotroph 
j = 6; %=>feed_rate = 7500 => Saturated effect;
sel2 = [4:8];
plot(auxfeed.Sin(sel2), reshape(mean(auxfeed.growth_rate(:,j,sel2),1),length(auxfeed.Sin(sel2)),[]),'--o','Displayname','Auxotroph','linewidth',2,'color',[.4 .4 .4])

% plot p = 3 equivalent to feed = 0; with lines
%j = 2;	% feed_rate = 312
j = 7;	% feed_rate = 12500

% Plot baseline with no feedback
k3_var = 1;
plot(Sin(2:end), reshape(mean(growth_rate(:,j,n_ext,2:end,k3_var),1),length(Sin(2:end)),[]),'--o','Displayname',strjoin({'% FdBk =',num2str(k3_var_range(k3_var))},' '),'linewidth',2,'color','k')


ct = 1;
%for k3_var = 1:length(k3_var_range);
for k3_var = [4:7];	% selected k3_var values
	plot(Sin(2:end), reshape(mean(growth_rate(:,j,n_ext,2:end,k3_var),1),length(Sin(2:end)),[]),'-o','Displayname',strjoin({'% FdBk =',num2str(k3_var_range(k3_var))},' '),'linewidth',2,'color',lin_colors(ct,:))
	ct = ct + 1;
end
legend('location','southeast','NumColumns',2);
%title({'Variation of growth rate with Substrate flux','and limiting metabolite feed rate',''});
ylabel('Mean growth rate (per hour)')
xlabel('Substrate flux available per cascade (per sec)')

f.Position(3:4) = f.Position(3:4)*0.75;
f.Children(2).XLabel.FontSize = 11;
f.Children(2).YLabel.FontSize = 11;

%name = 'Rev1_extFeed_GrowVsSin_feedback_lowFeed_312';
name = 'Rev1_extFeed_GrowVsSin_feedback_highFeed_12500';
set(f, 'Color', 'none');
export_fig(name, '-png', '-pdf','-r300', '-transparent', '-painters')




% Figure growth rate vs k3_var

fg = figure; hold on;
% plot p = 3 equivalent to feed = 0; with lines
%j = 2;
j = 7;
%j = 4;
for i = 1:length(Sin)
	plot(k3_var_range, reshape(mean(growth_rate(:,j,n_ext,i,:),1),length(k3_var_range),[]),'-o','Displayname',strjoin({'Sin =',num2str(Sin(i))},' '),'linewidth',2,'color',lin_colors(i,:))
end
legend('location','southeast','NumColumns',2);
%title({'Variation of growth rate with Substrate flux','and limiting metabolite feed rate',''});
ylabel('Mean growth rate (per hour)')
xlabel('% Feedback (Decrease in k3)')

%%===========================================
% Figure growth rate vs k3_var

plot_indices = [1 2 3 4 9];	% for feed_rate

%plot_indices = [1 2 3 4 5 9];	% for feed_rate

f = figure; hold on;

% plot p = 3 equivalent to feed = 0; with lines
i = 2;	% Sin = 5000
%i = 7;	% Sin = 13000

% Plot guide lines
%for j = 1:length(feed_rate(plot_indices))
%	plot(k3_var_range, repmat(mean(growth_rate(:,j,n_ext,i,1),1),1,length(k3_var_range)),'--','Displayname',strjoin({'Feed Rate =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',3,'color',[lin_colors(j,:) 0.5])
%end

for j = 1:length(feed_rate(plot_indices))
	plot(k3_var_range, reshape(mean(growth_rate(:,j,n_ext,i,:),1),length(k3_var_range),[]),'-o','Displayname',strjoin({'Feed Rate =',num2str(feed_rate(plot_indices(j)))},' '),'linewidth',2,'color',lin_colors(j,:))
end
legend('location','southwest','NumColumns',1);
%title({'Variation of growth rate with Substrate flux','and limiting metabolite feed rate',''});
ylabel('Mean growth rate (per hour)')
xlabel('% Feedback (Decrease in k3)')

f.Position(3:4) = f.Position(3:4)*0.5;
f.Children(2).XLabel.FontSize = 11;
f.Children(2).YLabel.FontSize = 11;

name = 'Rev1_extFeed_feedback_lowSin_5k';
%name = 'Rev1_extFeed_feedback_highSin_13k';
set(f, 'Color', 'none');
export_fig(name, '-png', '-pdf','-r300', '-transparent', '-painters')


