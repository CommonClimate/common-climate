function [field_r,diagn,st_i,mn_i] = mann09_cfr(field,proxy,t,calib,options)
% Function [field_r,diagn,st_i,mn_i] =  mann09_cfr(field,proxy,t,calib,options)
%
%    Emulates the nested hybrid RegEM TTLS CFR method of Mann et al (2009)
%
%  inputs:  - field, 2D climate field (nt x ns)
%           - proxy, proxy matrix (nt x nr)
%           - t, time axis  (nt)
%           - calib, index of calibration (< nt)
%           - options, structure. in particular:
%              # method = regem regularization scheme
%               ('ttls','iridge',mridge',ittls,'ipcr') [default = 'ttls']
%
%              # cutoff = period split (for hybrid RegEM only), e.g. 10 years
%              # nest = cell array containing nest indices as vectors  [default = 1:end]
%              # npcs  = Number of principal components used for compression (default = 4)
%              # splice = splicing mode, 'chrono' or 'RE' [default = 'chrono']
%
%  outputs: - field_r, reconstructed field
%           - diagn, structure of diagnostic outputs
%           -  mn_i, st_i: means and variances over the instrumental interval
%
%  ASSUMPTIONS:
%    -  Unlike most geologists, this code assumes that time runs FORWARD.
%
%  Reference: Mann, M. E., Zhang, Z., Rutherford, S., Bradley, R. S., Hughes, M. K., Shindell,
%    D., Ammann, C., Faluvegi, G., and Ni, F. (2009). Global signatures and dynamical
%   origins of the little ice age and medieval climate anomaly. Science, 326(5957):1256--1260.
%
%  Julien Emile-Geay, Written, 09/12/12, USC, based on regem_cfr_2p3 and hybrid_nest_regem_index.m
%     updated 08/26/2013, J.E.G 
%  >> uses hepta_smooth.m for filtering. 
%  >> no longer assumes that inputs are standardized (does so internally) 
% =================================================================================

% Process options
% ===============
if nargin < 3 | isempty(options)
   fopts = [];
else
   fopts = fieldnames(options);
end

% assign regularizatin method
if strmatch('method', fopts)  % number of principal components retain
   regmethod = options.method;
else
   regmethod = 'ttls';
end

if strmatch('regpar', fopts)  % number of principal components retain
   regpar = options.regpar;
else
   regpar = []; % by default, choose internally
end

% hybrid or not
if strmatch('cutoff', fopts)  % number of principal components retain
   try
      f_lw = 1/options.cutoff;
   catch
      warning('No cutoff period provided for hybrid split, defaults to 20 years');
      f_lw = 0.05;
   end
end


if strmatch('nest', fopts)  % minimum fraction of available values in a proxy segment
   nest = options.nest;
else
   nest{1}=[1:length(proxy)];
end

if strmatch('splice', fopts)
   splice = options.splice;
else
   splice = 'chrono';
end

%  Principal Component Compression
if strmatch('npcs', fopts)  % number of principal components retain
   npcs = options.npcs;
else
   npcs = [];
   %note: M09 used one different value per nest... see Table S3
end
%   ====== end options processsing  =======

% Define Time Parameters
% ========================
year_i = min(t);  year_f = max(t);
nr   = size(proxy,2); % number of proxy records
[ni,ns] = size(field); % field dimensions
nt = length(t);

% Reverse time for both field and proxies
%target.tot = field;
prox.tot   = proxy;
inst.tot   = field;

[inst_n.tot,mn_i.tot,st_i.tot] = standardize(inst.tot);

% Frequency split
display('Smoothing field');
for j=1:ns
   %  [inst.lw(:,j),icb,ice,mse0] = lowpassmin(inst.tot(:,j),f_lw);
   inst_n.lw(:,j) = hepta_smooth(inst_n.tot(:,j),f_lw);
end

display('Smoothing proxies');
prox.lw = NaN(size(prox.tot));
for k = 1:nr
   Z = prox.tot(:,k);
   %   prox.lw(~isnan(Z),k)  = lowpassmin(Z(~isnan(Z)),f_lw);
   prox.lw(~isnan(Z),k)  = hepta_smooth(Z(~isnan(Z)),f_lw);
end
display('Done smoothing');
% assign high-freq structures
inst_n.hi = inst_n.tot - inst_n.lw;
prox.hi = prox.tot - prox.lw;
%
st_p.hi = nstd(prox.hi(calib,:)); %st_i.hi = std(inst_n.hi);
st_p.lw = nstd(prox.lw(calib,:)); %st_i.lw = std(inst_n.lw);
st_p.tot = nstd(prox.tot(calib,:)); %st_i.tot = std(inst_n.tot);
%
mn_p.hi = nmean(prox.hi(calib,:)); %mn_i.hi = mean(inst_n.hi);
mn_p.lw = nmean(prox.lw(calib,:)); %mn_i.lw = mean(inst_n.lw);
mn_p.tot = nmean(prox.tot(calib,:)); %mn_i.tot = mean(inst_n.tot);

% Singular Value Decomposition of Instrumental Data
[Uh,Sh,Vh]      = svd(center(inst_n.hi),0); % "economy" decomposition
instpcs.hi      = Uh;
insteofs.hi     = Vh;
instsingvals.hi = Sh;
[Ul,Sl,Vl]      = svd(center(inst_n.lw),0); % "economy" decomposition
instpcs.lw      = Ul;
insteofs.lw     = Vl;
instsingvals.lw = Sl;
% input matrices are instrumental PCs
inst_mat.hi     = instpcs.hi;
inst_mat.lw     = instpcs.lw;
%np = npcs;  % number of variables to reconstruct is npcs

%==========================================================
% Millennnium Reconstruction
%==========================================================
N_t = length(nest);

for it = 1:N_t  % loop over nests
   % Find available proxy predictors
   inest = ismember(t,nest{it}); nit = sum(inest);
   cnest = ismember(nest{it},t(calib));
   P = prox.hi(inest,:); availn{it} = find(~isnan(std(P))); nprox(it) = numel(availn{it});
   % Assigning proxy matrices over this network
   prox_hi  = (prox.hi(inest,availn{it}) - repmat(mn_p.hi(availn{it}), [nit 1]))./repmat(st_p.hi(availn{it}), [nit 1]);
   prox_lw  = (prox.lw(inest,availn{it}) - repmat(mn_p.lw(availn{it}), [nit 1]))./repmat(st_p.lw(availn{it}), [nit 1]);
   prox_tot  = (prox.tot(inest,availn{it}) - repmat(mn_p.tot(availn{it}), [nit 1]))./repmat(st_p.tot(availn{it}), [nit 1]);
   
   % Assign number of PCs if needed
   if ~isempty('npcs')
      proxy_c =  prox_tot(cnest,:);
      [U,S,V] = svd(proxy_c,0); % do the decomposition, economy style. U's are the PCs
      d = diag(S).^2;
      dp = d./sum(d); % that gives the pct for each pc
      npcs(it) = eigenselect(dp)  % sthg is fucked up with my version of eigenselect
      %       % use Wax & Kailath method instead
      %       r = [1:50]';
      %       [wk85,ne08] = pca_truncation_criteria(dp, nr, r, ni, 1);
      %       [~,imin]=min(wk85.aicc);
      %       npcs(it) = r(imin)
   end
   np = npcs(it);
   
   % subsample the PC matrix according to this number
   inst_mat_hi = inst_mat.hi(:,1:np);
   inst_mat_lw = inst_mat.lw(:,1:np);
   
   switch regmethod
      case {'ttls'}  % if TTLS, need to determine truncation
         % with Mann 09 criterion
         if isempty(regpar)  % compute the regpar for hi&lw automatically
            % High-Frequency Truncation
            P = standardize([inst_mat_hi prox_hi(cnest,:)]);
            %P = [inst_mat.hi prox.hi(calib,:)];
            [U,S,V] = svd(P,0); % do the decomposition, economy style. U's are the PCs
            highfevalues = diag(S.^2);
            highfpctvar = diag(S).^2./sum(diag(S).^2); % that gives the pct for each pc
            mbest=eigenselect(highfpctvar);
            trunc.hi=mbest;
            %
            % Low-Frequency Truncation
            P = standardize([inst_mat_lw prox_lw(cnest,:)]);
            [U,S,V] = svd(P,0); % do the decomposition, economy style. U's are the PCs
            lowfevalues=diag(S.^2);
            lowfpctvar=diag(S).^2./sum(diag(S).^2); % that gives the pct for each pc
            cumlowfpctvar=cumsum(lowfpctvar);
            ttlslowfpar=length(find(cumlowfpctvar<=0.33))+1;
            trunc.lw=max(1,ttlslowfpar);
            %
            options.hi.regpar = trunc.hi;
            options.lw.regpar = trunc.lw;
         else  % specify the regpar as shown in the master code
            trunc.hi = options.hi.regpar;
            trunc.lw = options.lw.regpar;
         end
         display(trunc);
   end
   % Size
   [nit,pit]             =   size(prox_hi);
   %	Assemble climate field/proxy data matrices
   Xh_in				     = NaN(nit,np+pit);
   Xh_in(cnest,1:np)   = inst_mat_hi;  % Put in instrumental data
   Xh_in(:,np+1:np+pit) = prox_hi;  % the rest is (possibly incomplete) proxy data
   
   Xl_in			        = NaN(nit,np+pit);
   Xl_in(cnest,1:np)   = inst_mat_lw;  % Put in instrumental data
   Xl_in(:,np+1:np+pit) = prox_lw;  % the rest is (possibly incomplete) proxy data
   
   % Perform Hybrid RegEM (Parallel if nprocs > 1) . Use block-based version (if options.*.block = 1)
   disp('Handling high-freq recon');
   [Xh, Mh, Ch, Xerr_h, Wh] = regem_2p0(Xh_in,options.hi);
   disp('Handling low-freq recon');
   [Xl, Ml, Cl, Xerr_l, Wl] = regem_2p0(Xl_in,options.lw);
   
   %  SVD synthesis : field = U*S*V' = PCs*S*EOFs'; + restore mean and variance
   %recon{it}.hi  = Xh(:,1:np)*instsingvals.hi(1:np,:)*insteofs.hi'.*repmat(st_i.hi,[nit 1]);
   %recon{it}.hi  = recon{it}.hi + repmat(mn_i.hi,[nit 1]);
   recon.hi  = Xh(:,1:np)*instsingvals.hi(1:np,:)*insteofs.hi';
   %
   %recon{it}.lw  = Xl(:,1:np)*instsingvals.lw(1:np,:)*insteofs.lw'.*repmat(st_i.lw,[nit 1]);
   %recon{it}.lw  = recon{it}.lw + repmat(mn_i.lw,[nit 1]);
   recon.lw  = Xl(:,1:np)*instsingvals.lw(1:np,:)*insteofs.lw';
   %
   recon.tot = recon.hi + recon.lw;
   
   %field_r{it} = recon.tot;
   %  rescale and add mean back in 
   field_r{it} = recon.tot .* repmat(st_i.tot,[nit 1]) + repmat(mn_i.tot,[nit 1]);
   
   diagn{it}.err    = Xerr_l;
   diagn{it}.Ch     = Ch;
   diagn{it}.Cl     = Cl;
   diagn{it}.Mh     = Mh;
   diagn{it}.Ml     = Ml;
   diagn{it}.wl     = Wl;
   diagn{it}.wh     = Wh;
   diagn{it}.avail  = availn{it};
   diagn{it}.npcs   = npcs(it);
   diagn{it}.trunc  = trunc;
end


%field_r = flipdim(recon.tot,1); % convert to standard Earthian time
end








