function series_filtered = filterAnoms(series, RunMean, UpperLim)
% Inputs:
% RunMean = # of years to run filter to identify peaks outside "UpperLim" standard deviations

%% Standardize 
sigma = nanstd(series);
mu = nanmean(series);
series = (series - mu)./sigma;

seriesM = movmean(series,RunMean);
seriesInFill = movmean(series,5); %+/- 2.5 years infill
% index = abs(series-seriesM) > UpperLim;
index = abs(series-seriesM) > UpperLim*nanmean(abs(series-seriesM));
series(index) = seriesInFill(index);

% output
series_filtered = series.*sigma + mu;

end
