function [jet_lat_out, jet_speed_out, jet_width_out] = upscaleJetLat(jet_profile,lats,half_window,dlat,computewidth)
% function written by M. Osman (osmanm@mit.edu) May 2019 to upscale the
% jet latitude and speed to a finer resolution using an adaptive localized
% (i.e., defined by the centered half window) second-order polynomial fit
% INPUT
% jet_profile - an n x m matrix of n years specifiying the jet profile at m latitudes % e.g., NOAA20thCv2_avgLowTrop or WindSpeedAtl (iCESM)   
% lats - an m x 1 vector of the latitudes for which the jet profile is derived
% half_window - the number of nearby latitude cells to fit the polynomial to, defaults to 3
% dlat - the interpolation step for specifying the interpolated max lat/speed; defaults to 0.01;
% option - to compute the jet width following Barnes et al., 2013; defaults to false
warning('off','all');

if nargin < 3; half_window = 3; end
if nargin < 4; dlat = 0.01; end
if nargin < 5; computewidth = false; end
if size(jet_profile,2) ~= length(lats)
    warning('Warning! The number of latitude dimensions in ''jet_profile'' does not equal the latitude dimensions in ''lat''. Try again.');
    return;
end

jet_lat_in = nan(size(jet_profile,1));
jet_speed_in = nan(size(jet_profile,1));
for i = 1:length(jet_speed_in)
    [jet_speed_in(i), id] = nanmax(jet_profile(i,:));
    jet_lat_in(i) = lats(id);
end

jet_speed_out = nan(size(jet_profile,1),1);
jet_lat_out = nan(size(jet_speed_out));
    for j = 1:length(jet_lat_in)
        latRange = [round(jet_lat_in(j),2)-half_window, round(jet_lat_in(j),2)+half_window]; % p/m 5 degrees from naive lat
        % lats_interp(:,j) = [min(latRange):dlat:max(latRange)]';
        lats_interp = [min(latRange):dlat:max(latRange)]';
        indexer = lats >= min(latRange) & lats <= max(latRange);
        p = polyfit(lats(indexer),jet_profile(j,indexer)',2);
        % adaptive_jet_interp(:,j) = p(1).*lats_interp(:,j).^2 + p(2).*lats_interp(:,j).^1 + p(3).*lats_interp(:,j).^0;
        % [jet_speed_out(j,1), id] = nanmax(adaptive_jet_interp(:,j));
        % adaptive_jet_interp = p(1).*lats_interp(:,j).^2 + p(2).*lats_interp(:,j).^1 + p(3).*lats_interp(:,j).^0;        
        adaptive_jet_interp = p(1).*lats_interp.^2 + p(2).*lats_interp.^1 + p(3).*lats_interp.^0;        
        [jet_speed_out(j,1), id] = nanmax(adaptive_jet_interp); clearvars adaptive_jet_interp
        % jet_lat_out(j,1) = lats_interp(id,j);
        jet_lat_out(j,1) = lats_interp(id); clearvars lats_interp;
    end

if computewidth % calculate jet-width via; via Barnes et al., 2013, find width at values 1/2 the max speed, 
    jet_speed_interp_width = nan(size(jet_speed_out));
    jet_lat_interp_width = nan(size(jet_speed_out));
    jet_width_out = nan(size(jet_speed_out));
    lat_interp_width = [floor(max(lats)):-1*dlat:ceil(min(lats))]';
        for j = 1:length(jet_lat_in)
            jet_prof_width(:,j) = interp1(lats,jet_profile(j,:),lat_interp_width,'cubic');
            [jet_speed_interp_width(j,1), id] = nanmax(jet_prof_width(:,j));
            jet_lat_interp_width(j,1) = lat_interp_width(id);
        end
        for j = 1:length(jet_lat_in)
            lat_interp_width_upper = lat_interp_width(lat_interp_width>jet_lat_interp_width(j));
            lat_interp_width_lower = lat_interp_width(lat_interp_width<jet_lat_interp_width(j));
            [~, upper_jet_width_ind] = min(abs(jet_speed_interp_width(j)./2 - jet_prof_width(lat_interp_width>jet_lat_interp_width(j),j)));
            [~, lower_jet_width_ind] = min(abs(jet_speed_interp_width(j)./2 - jet_prof_width(lat_interp_width<jet_lat_interp_width(j),j)));
            jet_width_out(j,1) = abs(lat_interp_width_upper(upper_jet_width_ind) - lat_interp_width_lower(lower_jet_width_ind));
            clearvars lat_interp_width_upper lat_interp_width_lower
        end
end
    
if ~computewidth
    jet_width_out = [];
end
warning('on','all');

end