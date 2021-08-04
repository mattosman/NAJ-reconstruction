function [CVstats, CVstats_mu, CVstats_std, CVstats_med, CVstats_max] = KfoldCrossValJetInd(B,X,Ymu,Ysigma,JetObsLat,JetObsSpeed,JetObsWidth,lat,indexCal,indexVal)
% indices must be either "Lat", "Speed", or "Width"
% size of lat (l x 1) must be equal to size(B,2)

% this code performs K-fold Cross Validation for CCA-based regression
% See, e.g., Chapter 7 of Hastie, Tibshirani, Friedman, 2009: Elements 
% of Statistical Learning: Data Mining, Inference, and Prediction, or Faber
% and Rajko, Analytica Chimica Acta, 2007

% Written by M. Osman (osmanm@mit.edu) -- Apr 2020
% Computes the k-fold cross validation stats on just the jet-stream geometric indices (position, intensity, and width)

% INPUTS ------------------------------------------------------------------
% X = centered or standardized values of predictors, n samples by p predictors
% y = centered or standardized values of predictand, n samples by 1
% A = the number of modes to compute out to
% k = fold number, e.g., if 10 (default), regression is computed on 9/10 of the data,
%   and validation on remaining 1/10.  Is computed iteratively 10 times until
%   all data is used for validationf
% verbal = prints the integral time scale and number of cross validation tests 
%   into the command window; input either 'true' or 'false'; defaults to 'true'
% inflateVar = whether or not to inflate the variance in the cross validation
% OUTPUTS -----------------------------------------------------------------
% CVstats = a structure containing the following fields:
%   B: PLS slopes for all submodels 1:A (p x A x k)
%   f: PLS intercepts for all submodels 1:A (k x 1)
%   R2Calib: Coefficient of determination (% var explained) for y var's vs X var's over calibration interval (k x A)
%   r2Verif: Coefficient of determination (% var explained) for y var's vs X var's over verification interval (k x A)
%   RE: Reduction of Error, 1 - SSEfromobsy_verif / SSEfrommeany_calib (k x p)
%   CE: Coefficient of Error (k x p), 1 - SSEfromobsy_verif / SSEfrommeany_verif (k x A)
%   RMSEP: root mean square error of prediction (verification) (k x A)
% CVstats_mu and CVstats_std = two structures containing the mean and stdev of the above fields across all tests:
%   B: PLS slopes for all submodels 1:A (p x A)
%   f: PLS intercepts for all submodels 1:A (1xA)
%   R2Calib: Coefficient of determination (% var explained) for y var's vs X var's over calibration interval (1 x A)
%   r2Verif: Coefficient of determination (% var explained) for y var's vs X var's over verification interval (1 x A)
%   RE: Reduction of Error, 1 - SSEfromobsy_verif / SSEfrommeany_calib (1 x p)
%   CE: Coefficient of Error (k x p), 1 - SSEfromobsy_verif / SSEfrommeany_verif (1 x A)
%   RMSEP: root mean square error of prediction (verification) (1 x A)
% h0 = the length of the validation blocks (years)
% t0 = the integral timescale of y;
% =========================================================================

% run loop over all indices
for i = 1:size(indexVal,2)
    [CVstats.R2Calib(i,:), CVstats.r2Verif(i,:), CVstats.RE(i,:), CVstats.CE(i,:), CVstats.RMSEP(i,:)] =...
        CrossVal_R2RECE_jet(squeeze(B(:,:,i)),X,Ymu,Ysigma,JetObsLat,JetObsSpeed,JetObsWidth,indexVal(:,i),indexCal(:,i),lat); % each row represents a new cross-validatio (k), each column a different mode (A)
end

% Compute Mean and Stdev across all folds
    CVstats_mu.R2Calib  = nanmean(CVstats.R2Calib,1);  % R2Calib
    CVstats_med.R2Calib = nanmedian(CVstats.R2Calib,1);  % R2Calib
    CVstats_max.R2Calib = nanmax(CVstats.R2Calib,[],1);  % R2Calib
    CVstats_std.R2Calib = nanstd(CVstats.R2Calib,1,1); 
    CVstats_mu.r2Verif  = nanmean(CVstats.r2Verif,1);  % r2verif
    CVstats_med.r2Verif = nanmedian(CVstats.r2Verif,1);  % r2verif
    CVstats_max.r2Verif = nanmax(CVstats.r2Verif,[],1);  % r2verif
    CVstats_std.r2Verif = nanstd(CVstats.r2Verif,1,1); 
    CVstats_mu.CE       = nanmean(CVstats.CE,1);  % CE
    CVstats_med.CE      = nanmedian(CVstats.CE,1);  % CE
    CVstats_max.CE      = nanmax(CVstats.CE,[],1);  % CE
    CVstats_std.CE      = nanstd(CVstats.CE,1,1); 
    CVstats_mu.RE       = nanmean(CVstats.RE,1); % RE
    CVstats_med.RE      = nanmedian(CVstats.RE,1); % RE
    CVstats_max.RE      = nanmax(CVstats.RE,[],1); % RE
    CVstats_std.RE      = nanstd(CVstats.RE,1,1); 
    CVstats_mu.RMSEP    = nanmean(CVstats.RMSEP,1); % RMSEP
    CVstats_med.RMSEP   = nanmedian(CVstats.RMSEP,1); % RMSEP
    CVstats_max.RMSEP   = nanmax(CVstats.RMSEP,[],1); % RMSEP
    CVstats_std.RMSEP   = nanstd(CVstats.RMSEP,1,1); 

end

% =========================================================================
function [R2Calib,r2Verif,RE,CE,RMSEP] = CrossVal_R2RECE_jet(B,X,Ymu,Ysigma,JetObsLat,JetObsSpeed,JetObsWidth,indexVal,indexCal,lat)
% Written by M.Osman, Nov 2019 (osmanm@mit.edu)
% =========================================================================
    JetProfRec = X*B' .* Ysigma + Ymu ;
    [JetRecLat, JetRecSpeed, JetRecWidth] = upscaleJetLat(JetProfRec,lat,5,0.01,true);    
    % jet lat
    RMSEP(1)   = sqrt(nanmean((JetObsLat(indexVal)-JetRecLat(indexVal)).^2)); % root mean square error of prediction in verification period
    SSEP(1)    = sum((JetObsLat(indexVal)-JetRecLat(indexVal)).^2); % sum squared error of prediction in verification period
    Vc(1)      = sum((JetObsLat(indexVal)-mean(JetObsLat(indexCal))).^2); % sum squared error of calibration mean climatology
    Vv(1)      = sum((JetObsLat(indexVal)-mean(JetObsLat(indexVal))).^2); % sum squared error of verification mean climatology
    RE(1)      = 1-SSEP/Vc; % Reduction of Error statistic
    CE(1)      = 1-SSEP/Vv; % Cofficient of Error statistic 
    R2Calib(1) = corr(JetObsLat(indexCal),JetRecLat(indexCal))^2; % calibration squared pearson corr (observed vs. predictand)
    r2Verif(1) = corr(JetObsLat(indexVal),JetRecLat(indexVal))^2;   % verification squared pearson corr (observed vs. predictand)
    % jet speed
    RMSEP(2)   = sqrt(nanmean((JetObsSpeed(indexVal)-JetRecSpeed(indexVal)).^2)); % root mean square error of prediction in verification period
    SSEP(2)    = sum((JetObsSpeed(indexVal)-JetRecSpeed(indexVal)).^2); % sum squared error of prediction in verification period
    Vc(2)      = sum((JetObsSpeed(indexVal)-mean(JetObsSpeed(indexCal))).^2); % sum squared error of calibration mean climatology
    Vv(2)      = sum((JetObsSpeed(indexVal)-mean(JetObsSpeed(indexVal))).^2); % sum squared error of verification mean climatology
    RE(2)      = 1-SSEP/Vc; % Reduction of Error statistic
    CE(2)      = 1-SSEP/Vv; % Cofficient of Error statistic 
    R2Calib(2) = corr(JetObsSpeed(indexCal),JetRecSpeed(indexCal))^2; % calibration squared pearson corr (observed vs. predictand)
    r2Verif(2) = corr(JetObsSpeed(indexVal),JetRecSpeed(indexVal))^2;   % verification squared pearson corr (observed vs. predictand)
    % jet width
    RMSEP(3)   = sqrt(nanmean((JetObsWidth(indexVal)-JetRecWidth(indexVal)).^2)); % root mean square error of prediction in verification period
    SSEP(3)    = sum((JetObsWidth(indexVal)-JetRecWidth(indexVal)).^2); % sum squared error of prediction in verification period
    Vc(3)      = sum((JetObsWidth(indexVal)-mean(JetObsWidth(indexCal))).^2); % sum squared error of calibration mean climatology
    Vv(3)      = sum((JetObsWidth(indexVal)-mean(JetObsWidth(indexVal))).^2); % sum squared error of verification mean climatology
    RE(3)      = 1-SSEP/Vc; % Reduction of Error statistic
    CE(3)      = 1-SSEP/Vv; % Cofficient of Error statistic 
    R2Calib(3) = corr(JetObsWidth(indexCal),JetRecWidth(indexCal))^2; % calibration squared pearson corr (observed vs. predictand)
    r2Verif(3) = corr(JetObsWidth(indexVal),JetRecWidth(indexVal))^2;   % verification squared pearson corr (observed vs. predictand)    
end

% =========================================================================
function T0 = integtimsc(y)
% Source: Trenberth, 1984, Monthly Weather Review 112, p.2359
% Modified from Kinnard et al., 2011 by M. Osman, Jan 2019 (osmanm@mit.edu)
% INPUTS
%   y: a time series vector
% OUTPUTS
%   T0: the integral time scale of y, or the time between independent observations of y
    
% Initialize
    pk = nan(size(y)); % autocorrelation function of AR process
    pk(1,:) = 1; 
    N = size(pk,1);

% Calculate p(k) up to 1/4N lag
    pk = xcorr(detrend(y,0),round(length(y)/4),'biased');
    a = find(pk==max(pk));
    pk = pk(a:end)./pk(a);  % keep only one side of the correlogram, and scale p(0) to 1

    N = length(pk)-1;
% Integral time scale, Equation 2.8 Trentberth 1984    
    T0 = 1 + 2*sum((1-(1:N)'./N).*pk(2:end));   
end