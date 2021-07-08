%% Plot LOG GEV Graph with dataset points

figure;
% Plot the log of data as histogram
his = histogram(log(div_durs_exp/60),'Normalization','pdf');

% Plot the histogram data as scatter in original linear X Scale
%tst = figure; scatter(exp(his.BinEdges(2:end)),his.Values)

% histogram version
%tst = figure; histogram('BinEdges',exp(his.BinEdges),'BinCounts',his.Values,'FaceAlpha',0.5,'EdgeColor', 'none')
tst = figure; histogram('BinEdges',exp(his.BinEdges),'BinCounts',his.Values,'FaceAlpha',0.5,'EdgeColor', lines(1))

% Convert the Y Scale to log to show fit
tst.Children.YScale = 'log';
yl = ylim; xl = xlim

hold on;

% Fit the log data histogram to GEV, to obtain the log GEV fit.
pd = fitdist(log(div_durs_exp/60)','GeneralizedExtremeValue');

% Obtain the y values for the data for given x range
%x = 3.1:0.01:4;
x = log(xl(1)):0.01:log(xl(2));
y = pdf(pd, x);

% Plot the data on the graph
plot(exp(x),y,'linewidth',2);

ylim(yl)



