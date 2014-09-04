function Xb = block_bootstrap(X,len,Nb);
% function Xb = block_bootstrap(X,len,Nb);
%    Performs block-bootstrap resampling on timeseries X
%    X = original tiemseries
%    len = block length
%    Nb = number of desired boostrap samples 
% ============================================
nt = length(X);
ns = floor(nt/len);
Xb = zeros(nt,Nb); replace = 1;
% sample with replacement Nb times 
for b = 1:Nb
   sample = randsample(ns,ns,replace);
   for j = 1:ns
      Xb((j-1)*len+1:j*len,b) = X((sample(j)-1)*len+1:(sample(j)-1)*len+len);
   end
end

end