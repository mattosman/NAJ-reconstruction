function [Yrec, CCA, pctVar, CrossVal, TimeInt, pctVarp, pctVarj] = ccrNested(X,Y,smoother,ManAutComp,crossValSig)
% North Atlantic Jet Stream reconstruction using GrIS ice core records via CCR
% Written by M. Osman (osmanm@mit.edu), May 2019,
% 
% INPUTS:
%   Y - a MATLAB structure containing the variables "year" (n x 1) and "data" (n x M - standardized to cal. int.!), 
%       mu (Mx1; the calibration interval mean), sigma (Mx1; the calibration interval std)
%       var (Mx1; the calibration interval variance before standardization)
%   X - a MATLAB structure containing the variables "year" (N x 1) and "data" (N x P) and N > n
%       ** note that X and y do not have to be standardized or centered - code will do this;
%   smoother - an even integer over which to smooth the data series (if data is pre-smoothed, then
%       input "1"); default = 1
%   ManAutComp - decides whether to do manual or automatic variable selection; 
%       input either 'Man' or 'Aut'; defaults to 'Aut'
    if nargin < 4; ManAutComp = 'Aut';   end
    TF.Man = strcmp(ManAutComp,'Man');
    TF.Aut = strcmp(ManAutComp,'Aut');
    if ischar(ManAutComp) && (TF.Man ~= 1) && (TF.Aut ~= 1)
        warning('Please input either ''Man'' or ''Aut'' for the ''ManAutComp'' input and try again.'); return;
    elseif ~ischar(ManAutComp) 
        warning('Please input either ''Man'' or ''Aut'' for the ''ManAutComp'' input and try again.'); return; end
%   vFast - a 'true' or 'false' input; if true, computes the reconstruction
%       without the slow post hoc monte carlo cross validation significance tests; 
%       defaults to false
    if nargin < 5 
        crossValSig = true; 
    elseif ~islogical(crossValSig)
        error('Input for ''vFast'' is not a logical value! Try again.'); return;
    end
% OUTPUTS:
%   yrec: a structure containing the reconstructed series ("rec") and time
%       ("year") and uncertainty ("unc"; +/- 1 RMSEPcv)
%   CCA: an (m x 1) structure containing the final CCA output and input for each model (m)
%   pctVar: an (m x 1) cell showing amount of x (top row) and y-variance (bot row) 
%       described by calibration of model m
%   CrossVal: an (m x 1) structure containing following cross-validation details:
%       CVstats_mu (structure containing means's from k-CV procedures on observed data), 
%       CVstats_std (structure containing std's from k-CV procedures on observed data), 
%       CVstats_mc (statistical significance of CVstats_mu using pseudo-random surrogate data), 
%       k (# of CV folds per test), h (length (years) of CV folds), t (integral timescale of y, 
%       years), ManorAut (whether model run was created manually or automatically) and num_comp 
%       (number of orthogonal covariance modes ultimately chosen for prediction)
%   TimeInt: a (m x 2) vector containing the oldest and most recent year of model m

% Dependencies:
%   a folder named ccr with the following scripts:
%     ccr.m
%     ccrFast.m
%     ebisuzaki.m
%     invprctile.m
%     upscaleJetLat.m
%     ccr_nested.m
%     plotCrossValCCA.m
%     KfoldCrossValCCA.m
%     KfoldCrossValCCA_MC.m
%     KfoldCrossValCCA_BS.m
%     KfoldCrossValCCA_VA.m
%     KfoldCrossValJetInd_MC.m
%     KfoldCrossValJetInd.m

%% Check predictor/predictand arguments
if nanmax(Y.year) ~= nanmax(X.year)
    warning('Most recent year of predictor/predictand data are not equal! Try again.')
    return
end

% align both series so descending in time;
if Y.year(end) > Y.year(1)
    Y.year = flipud(Y.year); Y.data = flipud(Y.data);
end
if X.year(end) > X.year(1)
    X.year = flipud(X.year); X.data = flipud(X.data);
end

% define calibration period as the full y-series
calib_period = Y.year;
calib_index  = X.year >= min(calib_period) & X.year <= max(calib_period);

% default the smoother to 1 (no smoother applied)
if nargin < 3; smoother = 1; end

%% Fix the calibration intervals

proxyAvail = double(~isnan(X.data));
[~, col_ind] = find(proxyAvail); % row and column index
mat_ind = find(proxyAvail(:));  % matrix index
proxyAvail(mat_ind) = col_ind; % assign each column it's identifier
model_num = nan(size(proxyAvail,1),1);
for i = 1:size(proxyAvail)
    if i == 1
        model_num(i) = 1;
    else
        if isequal(proxyAvail(i,:), proxyAvail(i-1,:))
            model_num(i) = model_num(i-1);
        else
            model_num(i) = model_num(i-1) + 1;
        end
    end
end
num_models = [1:max(model_num)]';

%% Low pass smoothing option
if mod(smoother,2) == 0 && smoother >= 4 % if even and >= 4, apply smoother
    for i = 1:size(X.data,2)
        indexer = ~isnan(X.data(:,i));
        X.data(indexer,i) = lowpass(X.data(indexer,i),1/smoother,1,1);
    end
    for i = 1:size(Y.data,2)
        indexer = ~isnan(Y.data(:,i));
        Y.data(indexer,i) = lowpass(Y.data(indexer,i),1/smoother,1,1);
    end    
else
    warning('Input value for ''smoother'' is either odd-numbered or equal to 1. No smoothing applied.')
end

%% Preallocate output
Yrec.year = nan(size(X.data,1),1);
Yrec.rec = nan(size(X.data,1),size(Y.data,2));   %  Reconstructed JET
Yrec.unc = nan(size(X.data,1),size(Y.data,2));   %  Average model prediction error = cross-validation RMSE

CCA = cell(max(model_num),1);
pctVar = cell(max(model_num),1);
CrossVal = cell(max(model_num),1);
pctVarp = cell(max(model_num),1);
TimeInt = nan(max(model_num),2);

%% Run SVD to specify number of modes to test in the Y-data

thresh = 95;
[~,Sj,~] = svd(Y.data,'econ');
pctVarj = diag(Sj).^2/nansum(diag(Sj).^2).*100; pctVarj(:,2) = cumsum(pctVarj(:,1));
dj = nansum(double(pctVarj(:,2) <= thresh)) + 1;
disp(['Number of predictand components retained under ',num2str(thresh),'% retention criterion: ',num2str(dj),' modes.']);

%% Calibrate data

m = 1; % current model; initiate as model #1
    curr_row = find(model_num == m); % this variable is only for defining current predictor columns
    curr_col = find(proxyAvail(curr_row(1),:)); % simply look at the columns in the first row
    num_pred = length(curr_col); 

while m <= max(model_num) && num_pred > 1 % make sure there are multiple predictors
    disp(['Calibrating model ',num2str(m),' (',num2str(num_pred),' predictors)'])

    % define X_curr (prediction) and X_calib (calibration) for the current model calibration
%     X_curr.data  = X.data(1:max(current_row), current_col); % most recent years, down to current model's oldest year
%     X_curr.year  = X.year(1:max(current_row));
    X_curr.data  = X.data(curr_row, curr_col); % most recent years, down to current model's oldest year
    X_curr.year  = X.year(curr_row);
    X_calib.data = X.data(calib_index, curr_col); % calib_index previously defined for convienence; same in every model
    X_calib.year = X.year(calib_index);
    % Y_calib.data = Y.data; 
    % Y_calib.year = Y.year;
    
    % Save GBI and d18O series calibration interval means/stdevs for later
    X_calib.mean = nanmean(X_calib.data,1);  X_calib.stdev = nanstd(X_calib.data,0,1);
    X_curr.mean  = X_calib.mean;             X_curr.stdev = X_calib.stdev;  % for convenience
    % Y_calib.mean = nanmean(Y_calib.data,1);  Y_calib.stdev = nanstd(Y_calib.data,0,1);
    % now equalize all X data relavative to calibration interval by standardizing - leave y in real units!
    X_calib.data = (X_calib.data - X_calib.mean) ./ X_calib.stdev;
    X_curr.data  = (X_curr.data - X_curr.mean) ./ X_curr.stdev; % reminder: this the calibration interval mean/stdev
    % Y_calib.data = (Y_calib.data - Y_calib.mean); % retain variance in y, but center!
        
    % Calculate number of modes to retain for proxies:
    thresh_prox = 75;
    [~,Sp,~] = svd(X_calib.data,'econ');
    pctVarp_curr = diag(Sp).^2/nansum(diag(Sp).^2).*100; pctVarp_curr(:,2) = cumsum(pctVarp_curr(:,1));
    dp = nansum(double(pctVarp_curr(:,2) <= thresh_prox)) + 1;
    %    dp = round(dp/2);
    dcca = nanmin([dp, dj]);
    disp(['Number of predictor components retained under ',num2str(thresh_prox),'% retention criterion for model-',num2str(m),': ',num2str(dp),' modes.']);
    
    % =================== test optimal model =========================
    k = 2;
    h = waitbar(0,['Testing all canonical combinations for model ',num2str(m),'. Please wait...']);
    n = 1;
    for i = 2:dp  % changed from 1:dp
        waitbar(i/ dp);
        for j = 2:dj % changed from 1:d
            for g = 2:min([i,j]) % changed from 1:min([i,j])
                [~, CVstats_mu] = KfoldCrossValCCA(X_calib.data,Y.data,Y.mu,Y.sigma,i,j,g,k,false,Y.var); % weight by the variance in each jet lat
                MODES(n,:) = [i,j,g];
                CE(n,1) = CVstats_mu.CE_mu;
                RE(n,1) = CVstats_mu.RE_mu;
                RMSEP(n,1) = CVstats_mu.RMSEP_mu;
                n = n+1;
            end
        end
    end
    close(h); clearvars CVstats_mu
    
    % ================== define suggested number components ===============
    %   Goal: make a suggestion, based on a decision-tree; but allow the user
    %   to ultimately make the decision by showing the distribution of, 
    %   RMSEP, RE, and CE values as a function of mode
    [~, CVind_RMSEP] = min(RMSEP);
    [~, CVind_RE] = max(RE);
    [~, CVind_CE] = max(CE);
    sugg_dp = round(1/3*MODES(CVind_RMSEP,1) + 1/3*MODES(CVind_RE,1) + 1/3*MODES(CVind_CE,1)); 
    sugg_dj = round(1/3*MODES(CVind_RMSEP,2) + 1/3*MODES(CVind_RE,2) + 1/3*MODES(CVind_CE,2)); 
    sugg_dcca = round(1/3*MODES(CVind_RMSEP,3) + 1/3*MODES(CVind_RE,3) + 1/3*MODES(CVind_CE,3)); 
   
    if TF.Man == 1
    % If doing manual selection: plot the results to visualize
    plotCrossValCCA
        % Manually accept the results
        answer = questdlg(['Suggested number of retained components are as follows: dp = ', ...
            num2str(sugg_dp),'; dj = ', num2str(sugg_dj),'; dcc = ', num2str(sugg_dcca),'. Accept?'],...
        ['Choose Number of CCA Components for Model number ',num2str(m),'.'], ...
        'Yes','No','Yes');
        % Handle response
        switch answer
            case 'Yes'
                dp_input = sugg_dp; 
                dj_input = sugg_dj;
                dcca_input = sugg_dcca;
            case 'No'
                title = ['Enter number of components to retain for model number ',num2str(m),'.'];
                % dp
                prompt = ['Enter a number of dp components (between ',num2str(1),' and ',num2str(dp),').'];
                answer = inputdlg(prompt,title,[1 40],{num2str(sugg_dp)});
                dp_input = str2num(answer{1});
                if dp_input < 1 || dp_input > dp
                    warner = warndlg(['Input value must be between ',num2str(1),' and ',num2str(dp),...
                        '. Setting to suggested: ',num2str(sugg_dp)]);
                    dp_input = sugg_dp;
                    pause(2); close(warner); clearvars warner;
                end
                % dp
                prompt = ['Enter a number of dp components (between ',num2str(1),' and ',num2str(dj),').'];
                answer = inputdlg(prompt,title,[1 40],{num2str(sugg_dj)});
                dj_input = str2num(answer{1});
                if dj_input < 1 || dj_input > dj
                    warner = warndlg(['Input value must be between ',num2str(1),' and ',num2str(dj),...
                        '. Setting to suggested: ',num2str(sugg_dj)]);
                    dj_input = sugg_dj;
                    pause(2); close(warner); clearvars warner;
                end
                % dp
                prompt = ['Enter a number of dp components (between ',num2str(1),' and ',num2str(dcca),').'];
                answer = inputdlg(prompt,title,[1 40],{num2str(sugg_dcca)});
                dcca_input = str2num(answer{1});
                if dcca_input < 1 || dcca_input > dcca
                    warner = warndlg(['Input value must be between ',num2str(1),' and ',num2str(dcca),...
                        '. Setting to suggested: ',num2str(sugg_dcca)]);
                    dp_input = sugg_dcca;
                    pause(2); close(warner); clearvars warner;
                end                
        end; clearvars answer definput title
    close(CVfig1); close(CVfig2); close(CVfig3); % close figures until next loop;
    elseif TF.Aut == 1        
    	dp_input = sugg_dp; 
        dj_input = sugg_dj;
        dcca_input = sugg_dcca;
    end
    clearvars RE CE RMSEP MODES
    
    % ================= Conduct CCA regression analysis ===================
    disp(['Calibrating model ',num2str(m),': dp = ', num2str(dp_input),'; dj = ', num2str(dj_input),'; dcca = ', num2str(dcca_input),'.'])
    % bootstrap estimation
    BSiter = 300;
    U = nan(size(X_calib.data,2),size(X_calib.data,2),BSiter);
    S = nan(size(X_calib.data,2),size(X_calib.data,2),BSiter);
    V = nan(size(X_calib.data,1),size(X_calib.data,2),BSiter);
    [U(:,:,1),S(:,:,1),V(:,:,1)] = svd(X_calib.data','econ');
    for i = 2:BSiter
        ind = sort(randi([1, size(X_calib.data,2)], 1, size(X_calib.data,2))');
        [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(X_calib.data','econ');
        if i > 1
            [~,U(:,:,i)] = procrustes(U(:,:,1),U(:,:,i));
            [~,S(:,:,i)] = procrustes(S(:,:,1),S(:,:,i));
            [~,V(:,:,i)] = procrustes(V(:,:,1),V(:,:,i));
        end
    end 
    Up = nanmedian(U,3); Sp = nanmedian(S,3); Vp = nanmedian(V,3); 
    [Ut, St, Vt] = svd(Y.data','econ');
    [B] = ccrFast(Up,Sp,Vp,Ut,St,Vt,dp_input,dj_input,dcca_input); % Apply PLS on calibration block; 
    
    % ========================= CROSS-VALIDATION ==========================
    
    k = 2; % make sure this aligns with above
    [CVstats, CVstats_mu, CVstats_std, h, t, Bsave, indexVal,indexCal] = KfoldCrossValCCA_BS(X_calib.data,Y.data,Y.mu,Y.sigma,dp_input,dj_input,dcca_input,k,false,Y.var); % weight by the variance in each jet lat
    
    % upscale position, intensity, width
    [JetLat_obs, JetSpeed_obs, JetWidth_obs] = upscaleJetLat(Y.data .* Y.sigma + Y.mu,Y.lat,5,0.01,true);
    
    % run CE calculations for Jet Stream Indices (lat, speed, width)
    [~, CVstats_mu_JetInd, ~, CVstats_med_JetInd, CVstats_max_JetInd] = ...
        KfoldCrossValJetInd(Bsave,X_calib.data,Y.mu,Y.sigma,JetLat_obs,JetSpeed_obs,JetWidth_obs,Y.lat,indexCal,indexVal);    
    
    % reassign lat variables
    CVstats_mu_jetlat.R2Calib = CVstats_mu_JetInd.R2Calib(1);    CVstats_mu_jetlat.r2Verif = CVstats_mu_JetInd.r2Verif(1);     CVstats_mu_jetlat.CE = CVstats_mu_JetInd.CE(1);    CVstats_mu_jetlat.RE = CVstats_mu_JetInd.RE(1);    CVstats_mu_jetlat.RMSEP = CVstats_mu_JetInd.RMSEP(1);    
    CVstats_med_jetlat.R2Calib = CVstats_med_JetInd.R2Calib(1);  CVstats_med_jetlat.r2Verif = CVstats_med_JetInd.r2Verif(1);   CVstats_med_jetlat.CE = CVstats_med_JetInd.CE(1);  CVstats_med_jetlat.RE = CVstats_med_JetInd.RE(1);  CVstats_med_jetlat.RMSEP = CVstats_med_JetInd.RMSEP(1);    
    CVstats_max_jetlat.R2Calib = CVstats_max_JetInd.R2Calib(1);  CVstats_max_jetlat.r2Verif = CVstats_max_JetInd.r2Verif(1);   CVstats_max_jetlat.CE = CVstats_max_JetInd.CE(1);  CVstats_max_jetlat.RE = CVstats_max_JetInd.RE(1);  CVstats_max_jetlat.RMSEP = CVstats_max_JetInd.RMSEP(1);    
    % reassign speed variables
    CVstats_mu_jetspeed.R2Calib = CVstats_mu_JetInd.R2Calib(2);    CVstats_mu_jetspeed.r2Verif = CVstats_mu_JetInd.r2Verif(2);     CVstats_mu_jetspeed.CE = CVstats_mu_JetInd.CE(2);    CVstats_mu_jetspeed.RE = CVstats_mu_JetInd.RE(2);    CVstats_mu_jetspeed.RMSEP = CVstats_mu_JetInd.RMSEP(2);    
    CVstats_med_jetspeed.R2Calib = CVstats_med_JetInd.R2Calib(2);  CVstats_med_jetspeed.r2Verif = CVstats_med_JetInd.r2Verif(2);   CVstats_med_jetspeed.CE = CVstats_med_JetInd.CE(2);  CVstats_med_jetspeed.RE = CVstats_med_JetInd.RE(2);  CVstats_med_jetspeed.RMSEP = CVstats_med_JetInd.RMSEP(2);    
    CVstats_max_jetspeed.R2Calib = CVstats_max_JetInd.R2Calib(2);  CVstats_max_jetspeed.r2Verif = CVstats_max_JetInd.r2Verif(2);   CVstats_max_jetspeed.CE = CVstats_max_JetInd.CE(2);  CVstats_max_jetspeed.RE = CVstats_max_JetInd.RE(2);  CVstats_max_jetspeed.RMSEP = CVstats_max_JetInd.RMSEP(2);    
    % reassign width variables
    CVstats_mu_jetwidth.R2Calib = CVstats_mu_JetInd.R2Calib(3);    CVstats_mu_jetwidth.r2Verif = CVstats_mu_JetInd.r2Verif(3);     CVstats_mu_jetwidth.CE = CVstats_mu_JetInd.CE(3);    CVstats_mu_jetwidth.RE = CVstats_mu_JetInd.RE(3);    CVstats_mu_jetwidth.RMSEP = CVstats_mu_JetInd.RMSEP(3);    
    CVstats_med_jetwidth.R2Calib = CVstats_med_JetInd.R2Calib(3);  CVstats_med_jetwidth.r2Verif = CVstats_med_JetInd.r2Verif(3);   CVstats_med_jetwidth.CE = CVstats_med_JetInd.CE(3);  CVstats_med_jetwidth.RE = CVstats_med_JetInd.RE(3);  CVstats_med_jetwidth.RMSEP = CVstats_med_JetInd.RMSEP(3);    
    CVstats_max_jetwidth.R2Calib = CVstats_max_JetInd.R2Calib(3);  CVstats_max_jetwidth.r2Verif = CVstats_max_JetInd.r2Verif(3);   CVstats_max_jetwidth.CE = CVstats_max_JetInd.CE(3);  CVstats_max_jetwidth.RE = CVstats_max_JetInd.RE(3);  CVstats_max_jetwidth.RMSEP = CVstats_max_JetInd.RMSEP(3);    
    
    % ================== Monte Carlo tests of CV stats ====================
    
    if crossValSig
    
    % jet profile
    cv_mciter = 1e2;
    [CVstats_mcsig, mc_vals] = KfoldCrossValCCA_MC(X_calib.data,Y.data,Y.mu,Y.sigma,dp_input,dj_input,dcca_input,k, CVstats_mu, cv_mciter, 90);
    
    % jet indices
    cv_mciter = 1e2;
    [CVstats_mcsig_JetInd, ~] = KfoldCrossValJetInd_MC(Bsave,X_calib.data,Y.mu,Y.sigma,JetLat_obs,JetSpeed_obs,JetWidth_obs,Y.lat,indexCal,indexVal, CVstats_mu_JetInd, cv_mciter, 90); 
    
    % reassign lat variables
    CVstats_mcsig_jetlat.R2Calib_ObsPerc = CVstats_mcsig_JetInd.R2Calib_ObsPerc(1);  CVstats_mcsig_jetlat.r2Verif_ObsPerc = CVstats_mcsig_JetInd.r2Verif_ObsPerc(1);   CVstats_mcsig_jetlat.CE_ObsPerc = CVstats_mcsig_JetInd.CE_ObsPerc(1);  CVstats_mcsig_jetlat.RE_ObsPerc = CVstats_mcsig_JetInd.RE_ObsPerc(1);    CVstats_mcsig_jetlat.RMSEP_ObsPerc = CVstats_mcsig_JetInd.RMSEP_ObsPerc(1);    
    CVstats_mcsig_jetlat.R2Calib_threshlev = CVstats_mcsig_JetInd.R2Calib_threshlev(1);  CVstats_mcsig_jetlat.r2Verif_threshlev = CVstats_mcsig_JetInd.r2Verif_threshlev(1);   CVstats_mcsig_jetlat.CE_threshlev = CVstats_mcsig_JetInd.CE_threshlev(1);  CVstats_mcsig_jetlat.RE_threshlev = CVstats_mcsig_JetInd.RE_threshlev(1);    CVstats_mcsig_jetlat.RMSEP_threshlev = CVstats_mcsig_JetInd.RMSEP_threshlev(1);    
    CVstats_mcsig_jetlat.R2Calib_exceed = CVstats_mcsig_JetInd.R2Calib_exceed(1);  CVstats_mcsig_jetlat.r2Verif_exceed = CVstats_mcsig_JetInd.r2Verif_exceed(1);   CVstats_mcsig_jetlat.CE_exceed = CVstats_mcsig_JetInd.CE_exceed(1);  CVstats_mcsig_jetlat.RE_exceed = CVstats_mcsig_JetInd.RE_exceed(1);    CVstats_mcsig_jetlat.RMSEP_exceed = CVstats_mcsig_JetInd.RMSEP_exceed(1);    
    % reassign speed variables
    CVstats_mcsig_jetspeed.R2Calib_ObsPerc = CVstats_mcsig_JetInd.R2Calib_ObsPerc(2);  CVstats_mcsig_jetspeed.r2Verif_ObsPerc = CVstats_mcsig_JetInd.r2Verif_ObsPerc(2);   CVstats_mcsig_jetspeed.CE_ObsPerc = CVstats_mcsig_JetInd.CE_ObsPerc(2);  CVstats_mcsig_jetspeed.RE_ObsPerc = CVstats_mcsig_JetInd.RE_ObsPerc(2);    CVstats_mcsig_jetspeed.RMSEP_ObsPerc = CVstats_mcsig_JetInd.RMSEP_ObsPerc(2);    
    CVstats_mcsig_jetspeed.R2Calib_threshlev = CVstats_mcsig_JetInd.R2Calib_threshlev(2);  CVstats_mcsig_jetspeed.r2Verif_threshlev = CVstats_mcsig_JetInd.r2Verif_threshlev(2);   CVstats_mcsig_jetspeed.CE_threshlev = CVstats_mcsig_JetInd.CE_threshlev(2);  CVstats_mcsig_jetspeed.RE_threshlev = CVstats_mcsig_JetInd.RE_threshlev(2);    CVstats_mcsig_jetspeed.RMSEP_threshlev = CVstats_mcsig_JetInd.RMSEP_threshlev(2);    
    CVstats_mcsig_jetspeed.R2Calib_exceed = CVstats_mcsig_JetInd.R2Calib_exceed(2);  CVstats_mcsig_jetspeed.r2Verif_exceed = CVstats_mcsig_JetInd.r2Verif_exceed(2);   CVstats_mcsig_jetspeed.CE_exceed = CVstats_mcsig_JetInd.CE_exceed(2);  CVstats_mcsig_jetspeed.RE_exceed = CVstats_mcsig_JetInd.RE_exceed(2);    CVstats_mcsig_jetspeed.RMSEP_exceed = CVstats_mcsig_JetInd.RMSEP_exceed(2);    
    % reassign width variables
    CVstats_mcsig_jetwidth.R2Calib_ObsPerc = CVstats_mcsig_JetInd.R2Calib_ObsPerc(3);  CVstats_mcsig_jetwidth.r2Verif_ObsPerc = CVstats_mcsig_JetInd.r2Verif_ObsPerc(3);   CVstats_mcsig_jetwidth.CE_ObsPerc = CVstats_mcsig_JetInd.CE_ObsPerc(3);  CVstats_mcsig_jetwidth.RE_ObsPerc = CVstats_mcsig_JetInd.RE_ObsPerc(3);    CVstats_mcsig_jetwidth.RMSEP_ObsPerc = CVstats_mcsig_JetInd.RMSEP_ObsPerc(3);    
    CVstats_mcsig_jetwidth.R2Calib_threshlev = CVstats_mcsig_JetInd.R2Calib_threshlev(3);  CVstats_mcsig_jetwidth.r2Verif_threshlev = CVstats_mcsig_JetInd.r2Verif_threshlev(3);   CVstats_mcsig_jetwidth.CE_threshlev = CVstats_mcsig_JetInd.CE_threshlev(3);  CVstats_mcsig_jetwidth.RE_threshlev = CVstats_mcsig_JetInd.RE_threshlev(3);    CVstats_mcsig_jetwidth.RMSEP_threshlev = CVstats_mcsig_JetInd.RMSEP_threshlev(3);    
    CVstats_mcsig_jetwidth.R2Calib_exceed = CVstats_mcsig_JetInd.R2Calib_exceed(3);  CVstats_mcsig_jetwidth.r2Verif_exceed = CVstats_mcsig_JetInd.r2Verif_exceed(3);   CVstats_mcsig_jetwidth.CE_exceed = CVstats_mcsig_JetInd.CE_exceed(3);  CVstats_mcsig_jetwidth.RE_exceed = CVstats_mcsig_JetInd.RE_exceed(3);    CVstats_mcsig_jetwidth.RMSEP_exceed = CVstats_mcsig_JetInd.RMSEP_exceed(3);    
        
    end
    
    % ================ Model prediction + MC Uncertainty ==================
    
    % finally, need to recompute the uncertainty (cross-val RMSEP) by this time inflating the variables within 'KfoldCrossValPLS1.m': 
    [CVstats_RMSEP, CVstats_mu_RMSEP, CVstats_std_RMSEP, ~, ~] = KfoldCrossValCCA_VA(X_calib.data,Y.data,Y.mu,Y.sigma,dp_input,dj_input,dcca_input,k,false,Y.var); 
    CVstats.RMSEP = CVstats_RMSEP.RMSEP;
    CVstats_mu.RMSEP = CVstats_mu_RMSEP.RMSEP;
    CVstats_std.RMSEP = CVstats_std_RMSEP.RMSEP; 
        clearvars CVstats_RMSEP CVstats_mu_RMSEP CVstats_std_RMSEP
    Yrec.year_obs = Y.year;
    Yrec.obs = Y.data .* Y.sigma + Y.mu ;
    Yrec.year(curr_row) = X_curr.year;
    % Inflate the variance!! Added 10/31/19
        Yrecon = X_curr.data*B' .* Y.sigma + Y.mu ; % add back bias offset
        Yrecon_mean = nanmean(Yrecon,1);
        Yrecon = (Yrecon - Yrecon_mean) .* std(Yrec.obs,1,1)./std(Yrecon,1,1)  + Yrecon_mean;   
    Yrec.rec(curr_row,:) = Yrecon;    clearvars Yrecon Yrecon_mean;
    Yrec.unc(curr_row,:) = ones(size(Yrec.rec(curr_row,:))) .* CVstats_mu.RMSEP(1,:); % RMSEP over calibration interval    
    
    % ====================== Record Model-m output ========================
    
    CCA{m,1}.B = B;
    CCA{m,1}.Xcal = X_calib.data;
    for i = 1:size(Y.data,2)
        corr2(i) = corr(Yrec.obs(:,i),Yrec.rec(1:length(Yrec.year_obs),i)).^2;
    end
    pctVar{m,1} = corr2; clearvars corr2
    CrossVal{m,1}.CVstats_mu = CVstats_mu;
    CrossVal{m,1}.CVstats_std = CVstats_std;
    CrossVal{m,1}.CVstats_mu_jetlat = CVstats_mu_jetlat;     CrossVal{m,1}.CVstats_med_jetlat = CVstats_med_jetlat;     CrossVal{m,1}.CVstats_max_jetlat = CVstats_max_jetlat;
    CrossVal{m,1}.CVstats_mu_jetspeed = CVstats_mu_jetspeed; CrossVal{m,1}.CVstats_med_jetspeed = CVstats_med_jetspeed; CrossVal{m,1}.CVstats_max_jetspeed = CVstats_max_jetspeed;
    CrossVal{m,1}.CVstats_mu_jetwidth = CVstats_mu_jetwidth; CrossVal{m,1}.CVstats_med_jetwidth = CVstats_med_jetwidth; CrossVal{m,1}.CVstats_max_jetwidth = CVstats_max_jetwidth;
    if crossValSig
    CrossVal{m,1}.mc_vals = mc_vals;
    CrossVal{m,1}.CVstats_mcsig = CVstats_mcsig;
    CrossVal{m,1}.CVstats_mcsig_jetlat = CVstats_mcsig_jetlat;
    CrossVal{m,1}.CVstats_mcsig_jetspeed = CVstats_mcsig_jetspeed;
    CrossVal{m,1}.CVstats_mcsig_jetwidth = CVstats_mcsig_jetwidth;
    end
    CrossVal{m,1}.k = size(CVstats.B,3);
    CrossVal{m,1}.t = t;
    CrossVal{m,1}.h = h;
    CrossVal{m,1}.dp_input = dp_input; clearvars dp_input
    CrossVal{m,1}.dj_input = dj_input; clearvars dj_input
    CrossVal{m,1}.dcca_input = dcca_input; clearvars dcca_input
    pctVarp{m,1} = pctVarp_curr; clearvars pctVarp_curr
    TimeInt(m,:) = [nanmin(X_curr.year), nanmax(X_curr.year)];    

    % ===================== Define the next model =========================
    m = m+1; 
        if m > max(model_num)
            disp(['Out of (useful) proxies - model-building completed! ',num2str(m-1),'/',num2str(max(model_num)),' models completed.'])        
            break; 
        else
        curr_row = find(model_num == m); % this variable is only for defining current predictor columns
        curr_col = find(proxyAvail(curr_row(1),:)); % simply look at the columns in the first row
        num_pred = length(curr_col); 
        end
    
end

end % End of function