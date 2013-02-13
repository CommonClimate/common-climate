function [field,field_p,pattern,RE,CE,R2,contrib] = proxy_regress(X,M,C,i_ind,p_ind,calib)
% function [field,field_p,pattern,RE,CE,R2,contrib] = proxy_regress(X,M,C,i_ind,p_ind,calib)
% 		Estimates climate field from proxy records using specified values
% 		for the mean 'M' and covariance matrix 'C', e.g. from a bootstrap sample
%
% INPUTS:
%  - X: incomplete data matrix (n x p)
%  - M: mean (1 x p)
%  - C: covariance matrix (1 x p)
%  - i_ind, p_ind : indices of instrumental and proxy columns
%  - calib : calibration period
%
% OUTPUTS:
%  - field: completed field (n x pi)
%  - field_p: completed field (n x pi x np), one for each pattern of missing value
%  - pattern (cell array, length np)
%  - RE, CE, R2, verification statistics, one for each pattern of missing value
%  - contrb: matrix of proxy contributions to each point.
%
%   CAUTION: assumes that time runs FORWARD
%
%  Julien Emile-Geay, USC, 10/29/2012
%    coded while Hurricane Sandy was battering the eastern seaboard
% ======================================================================

[n,p] = size(X);
% get indices of missing values and initialize matrix of imputed values
indmis       = find(isnan(X));
nmis         = length(indmis);
if nmis == 0
    warning('No missing value flags found.')
    return                                      % no missing values
end
[jmis,kmis]  = ind2sub([n, p], indmis);
Xmis         = sparse(jmis, kmis, NaN, n, p); % matrix of imputed values

% for each row of X, assemble the column indices of the available
% values and of the missing values
kavlr        = cell(n,1);
kmisr        = cell(n,1);
for j=1:n
    kavlr{j}   = find(~isnan(X(j,:)));
    kmisr{j}   = find(isnan(X(j,:)));
end

% search for unique patterns
avail = double(~isnan(X));
[avail_uniq,I,J] = unique(avail,'rows','first');
np = size(avail_uniq,1); % number of patterns
Xp = zeros(n,p,np);
p_p = numel(p_ind); p_i = numel(i_ind);
contrib = zeros(n,p_i,p_p);

% remove mean
X = X - repmat(M, n, 1);

for j=1:np             % cycle over patterns
    % Fill in the pattern cell array
    avail_m = avail - repmat(avail_uniq(j,:),n,1);
    pattern{j} = find(std(avail_m,0,2) == 0)'; 
    jp     = min(pattern{j});             % position of this pattern
    pm     = length(kmisr{jp});  % number of missing values in this pattern
    mp     = length(pattern{j});
    if pm > 0
        % evaluate regression coefficients
        [B,S] = ggm(C,kavlr{jp},kmisr{jp});
        
        % missing value estimates
        Xmis(pattern{j}, kmisr{jp})  = X(pattern{j}, kavlr{jp}) * B;
        
        % make instrumental prediction
        Xp(:,kmisr{jp},j) = X(:, kavlr{jp}) * B;
        
        if nargout > 6
            %disp('Compute individual proxy contributions')
            for s = i_ind
                contrib(pattern{j},s,kavlr{jp}-p_i)  = X(pattern{j}, kavlr{jp}) .* repmat(B(:,s)',[mp 1]);
            end
        end
        
    else % instrumental period: use previous pattern
        %  CAUTION: this should be coded more generally, e.g. for cases where missing values
        %  do exist over the instumental period, in both proxy and instrumental matrices
        jp1     = min(pattern{j-1});
        % missing value estimates
        Xmis(pattern{j}, kmisr{jp1})  = X(pattern{j}, kavlr{jp1}) * B;
        % make instrumental prediction
        Xp(:,kmisr{jp1},j) = X(:, kavlr{jp1}) * B;
    end
end

% update data matrix X
X(indmis)  = Xmis(indmis);

% add mean to centered data matrix
X  = X + repmat(M, n, 1);

% add mean to centered data matrix
Xp = Xp + repmat(M, [n 1 np]);

% define output
field   = X(:,i_ind);
field_p = Xp(:,i_ind,:);

% compute verification statistics
ncal = sum(calib>0);
RE = zeros(n,p_i); R2 = RE; CE = RE;
field_c = field(calib,:);

for j=1:np
    [REn,CEn,R2n]=verif_stats(field_c,sq(field_p(calib,:,j)),1:ncal,1:ncal);
    RE(pattern{j},:) = repmat(REn,length(pattern{j}),1);
    CE(pattern{j},:) = repmat(CEn,length(pattern{j}),1);
    R2(pattern{j},:) = repmat(R2n,length(pattern{j}),1);
end

end
