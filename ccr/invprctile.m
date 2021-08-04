function p = invprctile(dist,y)
% inverse percentile
% INPUTS
%   dist = distribution of values (column vector); cannot contain NaN's!
%   y = value (scalar)
% OUTPUT
%   p = percentile of y in dist
% =========================================================================
sz    = size(dist);
dist  = sort(dist,1); 
q     = [0 100*(0.5:(sz(1)-0.5))./sz(1) 100]';
distx = [dist(1,:); dist(1:sz(1),:); dist(sz(1),:)];
p     = interp1q(distx,q,y); 
end
