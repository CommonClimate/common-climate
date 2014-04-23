function [snr,rho,signif,index,d] = snr_selection(lon_t,lat_t,lon_p,lat_p,inst,prox,options)

% This function is designed to compute the maximum value of correlation
% coefficient within a search radius.
% 
% USC Climate Dynamics group, 2012
% ============================================
if nargin < 7 || isempty(options)
    fopts  = [];
else
    fopts  = fieldnames(options);
end

if max(strcmp('sig_method', fopts)) ~= 1
    error('Significance test method must be specified.');
else
    display(['Using the ', options.sig_method, ' way to estimate the statistical significance of correlation'])
end

if max(strcmp('snr_method', fopts)) ~= 1
    snr_method = 'local';
else
    snr_method = options.snr_method;
    if strcmp(snr_method,'max')
        if strcmp('Rmax',fopts) ~= 1
            Rmax = 6370/3;
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
    if strcmp('local',snr_method)
        [~,index(k)]       = min(d(:,k));
        T                  = inst(:,index(k));
        [rho(k),signif(k)] = corr_sig(prox(:,k),T,options);
        
        % Maximized SNR
    elseif strcmp('max',snr_method)
        % select points within Rmax (search radius)
        ind{k}             = find(d(:,k) <= Rmax);
        T                  = inst(:,ind{k});
        % Compute correlations with those temperature points
        nT                 = length(ind{k});
        display([num2str(nT), ' temperature neighbors for proxy ', num2str(k)])
        R                  = corr(prox(:,k),T);
        [rho(k),index(k)]  = max(abs(R));
        [r(k),signif(k)]   = corr_sig(prox(:,k),T(:,index(k)),options);
        if r(k) - rho(k) > 10e-14
            error('rho is miscalculated')
        end
    end
    
end

snr = abs(rho)./sqrt(1-abs(rho).^2);

return
end