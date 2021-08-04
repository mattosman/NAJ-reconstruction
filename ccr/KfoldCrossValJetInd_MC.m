function [CVstats_mcsig, mc_vals] = ...
    KfoldCrossValJetInd_MC(B,X,Ymu,Ysigma,JetObsLat,JetObsSpeed,JetObsWidth,lat,indexCal,indexVal, CVstats_mu, num_tests, mc_thresh)
% Monte Carlo tests of CV stats for NAJ geometric indices

% Written by M. Osman (osmanm@mit.edu) Jan 2019

% Requires the "CVstats_mu" struct output from KfoldCrossValJetInd.m to run, 
% as well as the matlab functions "invprctile.m", and "ebisuzaki.m" to run.

% INPUTS ------------------------------------------------------------------
%   X -- the same X input to create the CVstats_mu struct
%   Y -- the same Y input to create the CVstats_mu struct
%   k -- the same k input to create the CVstats_mu struct
%   CVstats_mu -- mean value output from "KfoldCrossValPLS2.m" which forms
%       the observed values we wish to test against
%   CVstats_std -- std value output from "KfoldCrossValPLS2.m" which forms
%       the observed values we wish to test against
%   num_tests -- number of MC tests; defaults to 100
%   mc_thresh -- a number greater than 0 and less than 100; denotes the 
%       threshold percentile of each empirical distribution; defaults to 90
%   num_comp -- the number of components to include (must be < A)
% OUTPUTS -----------------------------------------------------------------
%   the CVstats_mcsig struct contains three variables:
%   XXX_ObsPerc = the observed percentile of the obverved CVstat value
%       against the MC-derived null distribution
%   XXX_MCthresh = the different MC values at the prescribed "mc_thresh"
%       value; note that the RMSEP threshold value is "100 - mc_thresh" since
%       we wish to minimize this value
%   XXX_exceed = a logical index determined by whether the observed CVstat
%       value exceeds the MCthresh value
% mc_vals -- the MC values refined in CVstats_mcsig
% =========================================================================

if nargin < 5
    num_tests = 1e2;
end
if nargin < 6
    mc_thresh = 90;
end

% standardize series
% jet_rec = (jet_rec - nanmean(jet_rec))./nanstd(jet_rec);
% jet_obs = (jet_obs - nanmean(jet_obs))./nanstd(jet_obs);

% preallocate:
	mc_vals.R2Calib_mc = nan(num_tests,3); 
	mc_vals.r2Verif_mc = nan(num_tests,3); 
	mc_vals.RE_mc      = nan(num_tests,3); 
	mc_vals.CE_mc      = nan(num_tests,3); 
	mc_vals.RMSEP_mc   = nan(num_tests,3); 
	h = waitbar(0,'Conducting jet indices significance tests.  Please wait...');
	for j = 1:num_tests
        waitbar(j / num_tests);
        surr_JetObsLat = ebisuzaki(JetObsLat);
        surr_JetObsSpeed = ebisuzaki(JetObsSpeed);
        surr_JetObsWidth = ebisuzaki(JetObsWidth);
        % surrogate_rec = surrogate_rec*std_jet_rec ./ sqrt(mean(surrogate_rec.^2));
        [~, CVstats_mu_mc, ~, ~, ~] = KfoldCrossValJetInd(B,X,Ymu,Ysigma,surr_JetObsLat,surr_JetObsSpeed,surr_JetObsWidth,lat,indexCal,indexVal); % update feb2020
        
        % store only the corresponding "num_comp" row
        mc_vals.R2Calib_mc(j,:) = CVstats_mu_mc.R2Calib;
        mc_vals.r2Verif_mc(j,:) = CVstats_mu_mc.r2Verif;
        mc_vals.RE_mc(j,:)      = CVstats_mu_mc.RE;
        mc_vals.CE_mc(j,:)      = CVstats_mu_mc.CE;
        mc_vals.RMSEP_mc(j,:)   = CVstats_mu_mc.RMSEP;
	end
    close(h);
                
    CVstats_mcsig.mc_thresh = mc_thresh;
        
        % compute inverse percentiles of the distributions
        for j = 1:3
            CVstats_mcsig.R2Calib_ObsPerc(j) = invprctile(mc_vals.R2Calib_mc(j,:), CVstats_mu.R2Calib(j));
                if isnan(CVstats_mcsig.R2Calib_ObsPerc(j)); CVstats_mcsig.R2Calib_ObsPerc(j) = 100 - 1/num_tests; end
            CVstats_mcsig.r2Verif_ObsPerc(j) = invprctile(mc_vals.r2Verif_mc(j,:), CVstats_mu.r2Verif(j));
                if isnan(CVstats_mcsig.r2Verif_ObsPerc(j)); CVstats_mcsig.r2Verif_ObsPerc(j) = 100 - 1/num_tests; end
            CVstats_mcsig.RE_ObsPerc(j)      = invprctile(mc_vals.RE_mc(j,:), CVstats_mu.RE(j));
                if isnan(CVstats_mcsig.RE_ObsPerc(j)); CVstats_mcsig.RE_ObsPerc(j) = 100 - 1/num_tests; end
            CVstats_mcsig.CE_ObsPerc(j)      = invprctile(mc_vals.CE_mc(j,:), CVstats_mu.CE(j));
                if isnan(CVstats_mcsig.CE_ObsPerc(j)); CVstats_mcsig.CE_ObsPerc(j) = 100 - 1/num_tests; end
            CVstats_mcsig.RMSEP_ObsPerc(j)   = invprctile(mc_vals.RMSEP_mc(j,:), CVstats_mu.RMSEP(j));
                if isnan(CVstats_mcsig.RMSEP_ObsPerc(j)); CVstats_mcsig.RMSEP_ObsPerc(j) = 0 + 1/num_tests; end
        end

        % compute inverse percentiles of the distributions - choose threshold; 
        CVstats_mcsig.R2Calib_threshlev = prctile(mc_vals.R2Calib_mc, mc_thresh);
        CVstats_mcsig.r2Verif_threshlev = prctile(mc_vals.r2Verif_mc, mc_thresh);
        CVstats_mcsig.RE_threshlev      = prctile(mc_vals.RE_mc, mc_thresh);
        CVstats_mcsig.CE_threshlev      = prctile(mc_vals.CE_mc, mc_thresh);
        CVstats_mcsig.RMSEP_threshlev   = prctile(mc_vals.RMSEP_mc, 100-mc_thresh); % want to minimize RMSE!
        
        % does the inverse percentile above exceed the threshold?
        CVstats_mcsig.R2Calib_exceed = CVstats_mcsig.R2Calib_ObsPerc > mc_thresh;
        CVstats_mcsig.r2Verif_exceed = CVstats_mcsig.r2Verif_ObsPerc > mc_thresh;
        CVstats_mcsig.RE_exceed      = CVstats_mcsig.RE_ObsPerc > mc_thresh;
        CVstats_mcsig.CE_exceed      = CVstats_mcsig.CE_ObsPerc > mc_thresh;
        CVstats_mcsig.RMSEP_exceed   = CVstats_mcsig.RMSEP_ObsPerc < 100 - mc_thresh;
            
end