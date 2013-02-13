function [spec,eig_vec,PC,RC,RCp,signif]=hepta_ssa(X,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function hepta_ssa : performs Singular Spectrum Analysis on time series X
%       with the method of Vautard and Ghil, Phys. D. 1989
%
% syntax : [spec,eig_vec,PC,RC,RCp] = hepta_ssa(X,[M, K])
%   Inputs 
%   X :  line vector
%   M : window length. Default value is M=length(X)/10.
%   K : number of EOFs used for reconstruction (0 by default)
%
% output : 
%   spec : eigenvalue spectrum, in % variance
%   eig_vec : eigenvector matrix ("temporal EOFs")
%   PC      : matrix of principal components
%   RC      : matrix of RCs (N*M, K) (only if K>0)
%   RCp     : Reconstructed timeseries
%   
%   Hepta Technologies, 2004
%    last updated 03/14/2012 to include automated choice for K (AICC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear signif
% first make sure the vector size is 1*N
if size(X,1)~= 1
    X=X';
end
N = length(X);
[Xr,mu,sigma] = standardize(X); % center the series
% set default value for M
if nargin('hepta_ssa')== 0
   M = round(N/10);
   K = 0;
elseif nargin('hepta_ssa')== -1 % nargin < 0 if variable # of args
   M = varargin{1};
   K = 0;
elseif nargin('hepta_ssa')== -2
   M = varargin{1};
   if ischar(varargin{2}) & strcmpi(varargin{2},'mcssa');
      mcssa = true;
      MC = 1000;
   else
      mcssa = false;
      K = varargin{2};
      signif = [1:K];
   end
end
Np=N-M+1;
% compute autocorelation
[gam,lags]=xcorr(Xr,M-1,'unbiased');
% Fill in Covariance matrix
C = toeplitz(gam(M:2*M-1)); %take postive half of autocorrelation diagram, hence M to 2M-1
% solve eigenvalue problem
[eig_vec, eig_val] = eigd(C);
spec    = eig_val/sum(eig_val);

% determine significant eigenvalues
if mcssa
   %[w,A,C,SBC,FPE,th] = arfit(Xr,1,1); %fit AR(1) model
   [a,var] = ar1(Xr); s = sqrt(var);
   noise = zeros(N,MC);
   noise(1,:) = repmat(Xr(1),[1 MC]);
   for jt=2:N  %
      noise(jt,:) = a*noise(jt-1,:)+ s*randn(1,MC);
   end
   noise = standardize(noise);
   for m = 1:MC
      [Gn,ln]       = xcorr(noise(:,m),M-1,'unbiased');
      Cn            = toeplitz(Gn(M:2*M-1)); 
      Lambda_R(:,m) =  diag(eig_vec * Cn * eig_vec'); % noise "eigenvalues"  
   end
   q95 = quantile(Lambda_R,0.95,2);
   signif = find(eig_val > q95)'; % index of modes rising above the background
   display(['MCSSA modes retained: ', int2str(signif)]);
   fig('MCSSA'),clf
   v = [1:M]'; ligr = [ 0.7000    0.7000    0.7000];
   lmin = min(Lambda_R,[],2); lmax = max(Lambda_R,[],2);
   area_fill(v',lmin',lmax',ligr,ligr,0,0.3),hold on
   plot(v,eig_val,'kx',v,q95,'r-','linewidth',[2])
    % try Kaiser rule?
elseif K == 0
    trunc = [0: length(spec)-1];
    [MDL, NE08, AIC, AICC] = pca_truncation_criteria(eig_val, 1, trunc, N, 1);
    [~,imin] = min(real(AICC));
    K = trunc(imin);
    display(['AICC truncation choice, K = ', int2str(K)]);
    signif = [1:K];
end


% compute PCs
decal=zeros(Np,M);
for t=1:N-M+1
    decal(t,:)=Xr(t:M+t-1);
end
PC = decal*eig_vec; % the columns of this matrix are Ak(t), k=1 to M

% compute reconstructed timeseries if K > 0
if ~isempty(signif)
    RC=zeros(N,length(signif));
     % first M terms   
     for t=1:M-1
        Av=flipud(PC(1:t,signif));
        eig_vec_red=eig_vec(1:t,signif);
        RC(t,:)=1/t*sum(Av.*eig_vec_red,1);    
     end
     %middle of timeseries
     for t=M:Np
        Av=flipud(PC(t-M+1:t,signif));
        eig_vec_red=eig_vec(1:M,signif);
        RC(t,:)=1/M*sum(Av.*eig_vec_red,1);    
     end  
     % Last M terms
     for t=Np+1:N
        Av=flipud(PC(t-M+1:Np,signif));
        eig_vec_red=eig_vec(t-N+M:M,signif);
        RC(t,:)=1/(N-t+1)*sum(Av.*eig_vec_red,1);  
     end  
     %sum and restore the mean and variance
     RCp = sigma*sum(RC,2)+ mu;
 else 
     RC = []; RCp = [];
end

