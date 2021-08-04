function [CVstats_mcsig, mc_vals] = KfoldCrossValCCA_MC(X, Y, Ymu, Ysigma, dp, dt, dcca, k, CVstats_mu, num_tests, mc_thresh)
% Monte Carlo tests of CV stats for NAJ profile

% Written by M. Osman (osmanm@mit.edu) Jan 2019

% Requires the "CVstats_mu" struct output from KfoldCrossValCCA.m to run, 
% as well as the matlab functions "KfoldCrossValCCA.m", "invprctile.m", and 
% and "ebisuzaki.m" to run.

% INPUTS ------------------------------------------------------------------
%   X -- the same X input to create the CVstats_mu struct
%   Y -- the same Y input to create the CVstats_mu struct
%   A -- the same A input to create the CVstats_mu struct
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

if nargin < 11
    num_tests = 1e2;
end
if nargin < 12
    mc_thresh = 90;
end
std_X = nanstd(X,1,1);

% preallocate:
	mc_vals.R2Calib_mc = nan(num_tests,size(Y,2)); 
	mc_vals.r2Verif_mc = nan(num_tests,size(Y,2)); 
	mc_vals.RE_mc      = nan(num_tests,size(Y,2)); 
	mc_vals.CE_mc      = nan(num_tests,size(Y,2)); 
	mc_vals.RMSEP_mc   = nan(num_tests,size(Y,2)); 
	h = waitbar(0,'Conducting jet profile significance tests.  Please wait...');
	for j = 1:num_tests
        waitbar(j / num_tests);
        % surrogate_Y = Y_surrogate(Y,false);
        surrogate_X = nan(size(X));
        for i = 1:size(X,2)
            surrogate_X(:,i) = ebisuzaki(X(:,i));
            surrogate_X(:,i) = surrogate_X(:,i)*std_X(i) ./ sqrt(mean(surrogate_X(:,i).^2));
        end
        % [~, CVstats_mu_mc] = KfoldCrossValPLS2(surrogate_X,surrogate_Y,A,k,false);
        % [~, CVstats_mu_mc, ~] = KfoldCrossValCCA(surrogate_X,surrogate_Y,Ymu,Ysigma,dp,dt,dcca,k,false); % original version;
        [~, CVstats_mu_mc, ~] = KfoldCrossValCCA(surrogate_X,Y,Ymu,Ysigma,dp,dt,dcca,k,false); % update feb2020
        % store only the corresponding "num_comp" row
        mc_vals.R2Calib_mc(j,:) = CVstats_mu_mc.R2Calib;
        mc_vals.r2Verif_mc(j,:) = CVstats_mu_mc.r2Verif;
        mc_vals.RE_mc(j,:)      = CVstats_mu_mc.RE;
        mc_vals.CE_mc(j,:)      = CVstats_mu_mc.CE;
        mc_vals.RMSEP_mc(j,:)   = CVstats_mu_mc.RMSEP;
	end
    close(h);
                
    CVstats_mcsig.mc_thresh = mc_thresh;
    for i = 1:size(Y,2)
        
        % compute inverse percentiles of the distributions
        CVstats_mcsig.R2Calib_ObsPerc(1,i) = invprctile(mc_vals.R2Calib_mc(:,i), CVstats_mu.R2Calib(1,i));
            if isnan(CVstats_mcsig.R2Calib_ObsPerc(1,i)); CVstats_mcsig.R2Calib_ObsPerc(1,i) = 100 - 1/num_tests; end
        CVstats_mcsig.r2Verif_ObsPerc(1,i) = invprctile(mc_vals.r2Verif_mc(:,i), CVstats_mu.r2Verif(1,i));
            if isnan(CVstats_mcsig.r2Verif_ObsPerc(1,i)); CVstats_mcsig.r2Verif_ObsPerc(1,i) = 100 - 1/num_tests; end
        CVstats_mcsig.RE_ObsPerc(1,i)      = invprctile(mc_vals.RE_mc(:,i), CVstats_mu.RE(1,i));
            if isnan(CVstats_mcsig.RE_ObsPerc(1,i)); CVstats_mcsig.RE_ObsPerc(1,i) = 100 - 1/num_tests; end
        CVstats_mcsig.CE_ObsPerc(1,i)      = invprctile(mc_vals.CE_mc(:,i), CVstats_mu.CE(1,i));
            if isnan(CVstats_mcsig.CE_ObsPerc(1,i)); CVstats_mcsig.CE_ObsPerc(1,i) = 100 - 1/num_tests; end
        CVstats_mcsig.RMSEP_ObsPerc(1,i)   = invprctile(mc_vals.RMSEP_mc(:,i), CVstats_mu.RMSEP(1,i));
            if isnan(CVstats_mcsig.RMSEP_ObsPerc(1,i)); CVstats_mcsig.RMSEP_ObsPerc(1,i) = 0 + 1/num_tests; end

        % compute inverse percentiles of the distributions - choose threshold; 
        CVstats_mcsig.R2Calib_threshlev(1,i) = prctile(mc_vals.R2Calib_mc(:,i), mc_thresh);
        CVstats_mcsig.r2Verif_threshlev(1,i) = prctile(mc_vals.r2Verif_mc(:,i), mc_thresh);
        CVstats_mcsig.RE_threshlev(1,i)      = prctile(mc_vals.RE_mc(:,i), mc_thresh);
        CVstats_mcsig.CE_threshlev(1,i)      = prctile(mc_vals.CE_mc(:,i), mc_thresh);
        CVstats_mcsig.RMSEP_threshlev (1,i)  = prctile(mc_vals.RMSEP_mc(:,i), 100-mc_thresh); % want to minimize RMSE!
        
        % does the inverse percentile above exceed the threshold?
        CVstats_mcsig.R2Calib_exceed(1,i) = CVstats_mcsig.R2Calib_ObsPerc(1,i) > mc_thresh;
        CVstats_mcsig.r2Verif_exceed(1,i) = CVstats_mcsig.r2Verif_ObsPerc(1,i) > mc_thresh;
        CVstats_mcsig.RE_exceed(1,i)      = CVstats_mcsig.RE_ObsPerc(1,i) > mc_thresh;
        CVstats_mcsig.CE_exceed(1,i)      = CVstats_mcsig.CE_ObsPerc(1,i) > mc_thresh;
        CVstats_mcsig.RMSEP_exceed(1,i)   = CVstats_mcsig.RMSEP_ObsPerc(1,i) < 100 - mc_thresh;
        
    end
    
end