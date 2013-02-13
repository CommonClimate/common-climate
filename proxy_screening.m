function [rho,signif,index,d] = proxy_screening(lon_t,lat_t,lon_p,lat_p,inst,prox,options)
%  Screen for significant correlations between a proxy database and an instrumental field.  
%  The target temperature is either the closest gridpoint ('local') or maximum correlations 
%   within a specified search radius (default = 6370 km).
%   
%   The core of the problem is testing for correlations between timeseries that violate
%    the assumption of IID (indepedent and identically distributed) random variables.
%    The first "I" is the issue: most geophysical timeseries exhibit important serial
%    correlation, which violates indepdence over time. There are 3 ways to correct for this:
%
%     - ['ttest'] correct the classical T-test for the reduction of degrees of freedom induced by autocorrelation 
%     - ['isopersistent']: estimate the null distrubution by modeling the timeseries as AR(1) processes 
%    -  ['isospectral']: estimate the null distrubution by phase randomization
%
%    The function corr_sig.m performs this task, with options.
%
%    Note that all methods still require the "ID" assumption, specifically that of **Gaussianity**. 
%    Make sure your proxies are all reasonably close to a normal distribution before screening
%    for significant correlations (for instance, by applying boxcox_hepta.m)
%
%   INPUTS:
%	    - lon_t, lat_t (vectors): longitude and latitude of temperature gridpoints
%      - lon_p, lat_p (vectors): longitude and latitude of proxy locations
%      - inst, prox: instrumental and  
%   OUTPUTS:
%      - rho: correlation between each proxy and chosen temperature point
% 
%   N.B: make sure both matrices are correctly aligned (i.e. temperature and proxies 'see' the same time)
%         and that longitudes are in the [-180;+180] interval.
%
%  USC Climate Dynamics Lab, 2012.
%   based on snr_selection.m by Jianghao Wang.
% ==============================================

if nargin < 7 || isempty(options)
    fopts  = [];
else
    fopts  = fieldnames(options);
end

if max(strcmp('sig_method', fopts)) ~= 1
    display('Default significance test: ISOSPECTRAL')
else
    display(['Using the ', options.sig_method, ' way to estimate the statistical significance of correlation'])
end

if max(strcmp('screening_method', fopts)) ~= 1
    screening_method = 'local';
else
    screening_method = options.screening_method;
    if strcmp(screening_method,'max')
        if strcmp('Rmax',fopts) ~= 1
            Rmax = 6370;
        else
            Rmax = options.Rmax;
        end
    end
end

nr     = length(lon_p);ns = length(lon_t);
index  = zeros(nr,1); ind = cell(nr,1);
rho    = zeros(nr,1); r   = zeros(nr,1);
signif = zeros(nr,1);
d      = zeros(ns,nr);
c      = pi/180;

for k = 1:nr
    % Compute distance between each proxy point and HadCRUT3v grid point
    for i = 1:ns
        d(i,k) = greatCircleDistance(c*lat_t(i), c*lon_t(i), c*lat_p(k), c*lon_p(k));
    end
    
    % local SNR
    if strcmp('local',screening_method)
        [~,index(k)]       = min(d(:,k));
        T                  = inst(:,index(k));
        [rho(k),signif(k)] = corr_sig(prox(:,k),T,options);
        
        % Maximized SNR
    elseif strcmp('max',screening_method)
        % select points within Rmax (search radius)
        ind{k}             = find(d(:,k) <= Rmax);
        T                  = inst(:,ind{k});
        % Compute correlations with those temperature points
        nT                 = length(ind{k});
        display([num2str(nT), ' temperature neighbors for proxy ', num2str(k)])
        R                  = corr(prox(:,k),T);
        [rho(k),index(k)]  = max(abs(R));
        [r(k),signif(k)]   = corr_sig(prox(:,k),T(:,index(k)),options);
        %if abs(r(k) - rho(k)) > 10e-14
        %    error('rho is miscalculated')
        %end
    end
    
end

signif = logical(signif);
return
end
