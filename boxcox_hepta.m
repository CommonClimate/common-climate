function [z,lambda] = boxcox_hepta(y);

% function [z,lambda]=boxcox_hepta(y);
% Box-Cox transformation for a univariate series.
%
% input  -	y : data vector (positive values only - if not, shift upwards)
%
% output 
%		 -  lambda : maximum-likelhood estimate of lambda
%		 -  z      : Box-Cox transformed version of y, verifying:
%
%              | (y^lambda-1)/lambda  if lambda ~= 0 
%          z = | 
%              |  log(y)   otherwise 
%
%   Note: the numerical tolerance for lambda is set at 1e-5
%
%  Julien Emile-Geay, USC, June 2012
%  
% Reference: 
%   Box, G., and D. Cox (1964), An analysis of transformations, J. Roy. Stat. Soc., Ser. B
%    
% ===================================================================

y    = y(:); %ensure that y is a column vector and positive definite
ymin = nmin(y);
n = length(y);

% if ymin < 0
%    error('The Box-Cox transformation is only defined on variables taking positive values \n')
%    display('If your dataset crosses zero, simply apply the function to y - min(y)')
% end

%ensure that y is positive definite and free of NaNs
yn = y(~isnan(y))- ymin + eps;

% find optimal exponent
minopt = optimset('TolX', 1e-4, 'Display', 'off');
lambda  = fminbnd(@boxcox_llh, -2, 2, minopt, yn);

% apply power transform
zn = power_transform(lambda,yn);
z = NaN(n,1);
z(~isnan(y)) = zn; % assign the transform to all defined values
 
end

% ============= Box Cox log-likelihood =============
function llh = boxcox_llh(a,x)
   xt = power_transform(a,x);
   n = length(x); nh = n/2;
   v = var(xt);  
   % compute -loglik, so that its minimum coincides with maximum likelihood.
   llh = nh*log(v)-(a-1).*sum(log(x)) ;  
end

% ============= power transform =============
function xt = power_transform(a,x)
if abs(a) < 1e-5
   xt = log(x);
else
   xt = (x.^a-1)./a;
end
end