function [CVstats, CVstats_mu, CVstats_std, h0, t0, Bsave, indexVal,indexCal] = ...
KfoldCrossValCCA_BS(X,Y,Ymu,Ysigma,dp,dt,dcca,k,verbal,weights)
% this code performs K-fold Cross Validation for CCA-based regression
% See, e.g., Chapter 7 of Hastie, Tibshirani, Friedman, 2009: Elements 
% of Statistical Learning: Data Mining, Inference, and Prediction, or Faber
% and Rajko, Analytica Chimica Acta, 2007

% KfoldCrossValCCA version BootStrap (BS): 
% Adds a bootstrap estimation routine to KfoldCrossValCCA.m - updated Mar2020

% Written by M. Osman (osmanm@mit.edu) -- May 2019
% Inspired by moveblockCVpls.m presented by Kinnard et al., 2011 Nature

% INPUTS ------------------------------------------------------------------
% X = centered or standardized values of predictors, n samples by p predictors
% Y = centered or standardized values of predictand, n samples by m
% Ymu = mean values of the Yvariable w/ lat (mx1)
% Ysigma = standard devation values of the Yvariable w/ lat (mx1)
% dp = number of retained proxy modes
% dt = number of retained predictand modes
% dcca = number of retained cca modes
% k = fold number, e.g., if 10 (default), regression is computed on 9/10 of the data,
%   and validation on remaining 1/10.  Is computed iteratively 10 times until
%   all data is used for validationf
% verbal = prints the integral time scale and number of cross validation tests 
%   into the command window; input either 'true' or 'false'; defaults to 'true'
% weights for summing the calibration statistics; could use, e.g., the
%   latitudinal area, or the variance of each grid
% OUTPUTS -----------------------------------------------------------------
% CVstats = a structure containing the following fields:j
%   B: PLS slopes for all submodels 1:A (p x A x k)
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

if nargin < 8; k = 10; end % default value
if nargin < 9; verbal = true; end
if nargin < 10; weights = ones(1,size(Y,2)); end

% compute integral timescale
T0i = nan(size(Y,2),1); % preallocate
for i = 1:size(Y,2)
	T0i(i,1) = round(integtimsc(Y(:,i)));
    if T0i(i,1) == 0; T0i(i,1) = T0i(i,1)+1; end
end
T0 = round(nanmean(T0i)); % take mean of the integral timescale of all lat-vals
t0 = T0; % for output

h = floor(length(Y)/k); h0 = h;
% steps by integral timescale, so all points are independent
    indexVal = zeros(length(Y),1);
    i = 1;
    j = 1;
    while h <= size(indexVal,1)
        indexVal(i:h,j) = 1;
        j = j+1;
        i = i + T0; 
        h = h + T0;
    end
    indexVal = logical(indexVal);
    indexCal = zeros(size(indexVal));
    indexCal(~indexVal) = 1; indexCal = logical(indexCal);     
    if verbal
        disp(['Integral timescale: ',num2str(T0),' --- Number of Cross Validation tests: ',num2str(size(indexVal,2))]);
    end
    
BSiter = 300;
% run loop over all indices
for i = 1:size(indexVal,2)
    XCal = X(indexCal(:,i),:); % calibration X block
    YCal = Y(indexCal(:,i),:);   % calibration Y block  
    % bootstrap estimation
    U = nan(size(XCal,2),size(XCal,2),BSiter);
    S = nan(size(XCal,2),size(XCal,2),BSiter);
    V = nan(size(XCal,1),size(XCal,2),BSiter);
    [U(:,:,1),S(:,:,1),V(:,:,1)] = svd(XCal','econ');
    for j = 2:BSiter
        ind = sort(randi([1, size(XCal,2)], 1, size(XCal,2))');
        [U(:,:,j),S(:,:,j),V(:,:,j)] = svd(XCal','econ');
        if j > 1
            [~,U(:,:,j)] = procrustes(U(:,:,1),U(:,:,j));
            [~,S(:,:,j)] = procrustes(S(:,:,1),S(:,:,j));
            [~,V(:,:,j)] = procrustes(V(:,:,1),V(:,:,j));
        end
    end 
    Up = nanmedian(U,3); Sp = nanmedian(S,3); Vp = nanmedian(V,3); 
    [Ut, St, Vt] = svd(YCal','econ');
    [B] = ccrFast(Up,Sp,Vp,Ut,St,Vt,dp,dt,dcca); % Apply PLS on calibration block; 
    CVstats.B{i} = B; % Compute and store CV statistics
    [CVstats.R2Calib(i,:), CVstats.r2Verif(i,:), CVstats.RE(i,:), CVstats.CE(i,:), CVstats.RMSEP(i,:)] =... % (m x A x k), where k is the CV iteration
        CCA_CrossVal_R2RECE(X,Y,Ymu,Ysigma,B,indexVal(:,i),indexCal(:,i)); % each row represents a new cross-validatio (k), each column a different mode (A)
    Bsave(:,:,i) = B;
end

% Compute Mean and Stdev across all folds 
    Bdummy = nan(size(B,1),size(B,2),size(indexVal,2));
        for i = 1:size(indexVal,2) % loop through each CV, grab the jth layer
        	Bdummy(:,:,i) = CVstats.B{i};
        end
    CVstats_mu.B = nanmean(Bdummy,3);  % PLS slope standard error % ---- need to loop to save mean of each layer
    CVstats_std.B = nanstd(Bdummy,1,3); 
    CVstats_mu.R2Calib  = nanmean(CVstats.R2Calib,1);  % R2Calib
    CVstats_std.R2Calib = nanstd(CVstats.R2Calib,1,1); 
    CVstats_mu.r2Verif  = nanmean(CVstats.r2Verif,1);  % r2verif
    CVstats_std.r2Verif = nanstd(CVstats.r2Verif,1,1); 
    CVstats_mu.CE       = nanmean(CVstats.CE,1);  % CE
    CVstats_std.CE      = nanstd(CVstats.CE,1,1); 
    CVstats_mu.RE       = nanmean(CVstats.RE,1); % RE
    CVstats_std.RE      = nanstd(CVstats.RE,1,1); 
    CVstats_mu.RMSEP    = nanmean(CVstats.RMSEP,1); % RMSEP
    CVstats_std.RMSEP   = nanstd(CVstats.RMSEP,1,1); 
    
    % store the weighted mean across all latitudes
    CVstats_mu.R2Calib_mu  = nansum(weights .* CVstats_mu.R2Calib) ./ nansum(weights);  % R2Calib 
    CVstats_std.R2Calib_mu = nansum(weights .* CVstats_std.R2Calib) ./ nansum(weights);   
    CVstats_mu.r2Verif_mu  = nansum(weights .* CVstats_mu.r2Verif) ./ nansum(weights);   % r2verif
    CVstats_std.r2Verif_mu = nansum(weights .* CVstats_std.r2Verif) ./ nansum(weights);  
    CVstats_mu.CE_mu       = nansum(weights .* CVstats_mu.CE) ./ nansum(weights);    % CE
    CVstats_std.CE_mu      = nansum(weights .* CVstats_std.CE) ./ nansum(weights); 
    CVstats_mu.RE_mu       = nansum(weights .* CVstats_mu.RE) ./ nansum(weights);  % RE
    CVstats_std.RE_mu      = nansum(weights .* CVstats_std.RE) ./ nansum(weights); 
    CVstats_mu.RMSEP_mu    = nansum(weights .* CVstats_mu.RMSEP) ./ nansum(weights);  % RMSEP
    CVstats_std.RMSEP_mu   = nansum(weights .* CVstats_std.RMSEP) ./ nansum(weights);  
    
end

% =========================================================================
function [R2Cal,r2Val,RE,CE,RMSEP] = CCA_CrossVal_R2RECE(X,Y,Ymu,Ysigma,B,indexVal,indexCal)
% Written by M.Osman, Jan 2019 (osmanm@mit.edu)
% Inputs:
%   X -- matrix of standardized predictors (n x p)
%   Y -- vector of centered (or standardized, optional) predictand (n x m)
%   B -- PLS1 coefficients for the calibration period ([p+1] x m x A)
%   indexVer -- logical indexer denoting the validation period (n x 1)
%   indexCal -- logical indexer denoting the calibration period (n x 1)
% =========================================================================
    Ytrue  = Y .* Ysigma + Ymu;
    YValid = X(indexVal,:)*B' .* Ysigma + Ymu;     % regression prediction for verification period
    YCalib = X(indexCal,:)*B' .* Ysigma + Ymu;  
    RMSEP  = sqrt(nanmean((Ytrue(indexVal,:)-YValid).^2)); % root mean square error of prediction in verification period
    SSEP   = sum((Ytrue(indexVal,:)-YValid).^2); % sum squared error of prediction in verification period
    Vc     = sum((Ytrue(indexVal,:)-mean(Ytrue(indexCal,:))).^2); % sum squared error of calibration mean climatology
    Vv     = sum((Ytrue(indexVal,:)-mean(Ytrue(indexVal,:))).^2); % sum squared error of verification mean climatology
    RE     = 1-SSEP./Vc;  % Reduction of Error statistic
    CE     = 1-SSEP./Vv;  % Cofficient of Error statistic 
    for i = 1:size(Y,2)
        R2Cal(1,i) = corr(Ytrue(indexCal,i),YCalib(:,i)).^2; % calibration squared pearson corr (observed vs. predictand)
        r2Val(1,i) = corr(Ytrue(indexVal,i),YValid(:,i)).^2;   % verification squared pearson corr (observed vs. predictand)   
    end    

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