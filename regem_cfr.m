function [field_r,diagn] = regem_cfr(field,proxy,t,calib,options)
% Function [field_r,diagn] = regem_cfr(field,proxy,t,calib,options)
%  assess quality of RegEM reconstruction of surface temperature field in pseudoproxy context
%
%  inputs:  - field, 2D climate field (nt x ns)
%           - proxy, proxy matrix (nt x nr)
%           - t, time axis  (nt)
%           - calib, index of calibration (< nt)
%           - options, structure. in particular:
%              # method = regem regularization scheme
%               ('ttls','iridge',mridge',ittls,'ipcr')
%              # hybrid = boolean variable (true --> hybrid RegEM)
%              # cutoff = period split (for hybrid RegEM only)
%              # all others are standard RegEM options
%
%  outputs: - field_r, reconstructed field
%           - diagn, structure of diagnostic outputs
%
%  Written, 08/24/11, USC, based on cfr_pproxy_recon_regem_hybrid
%    main diffrences:
%     1) removed the pseudoproxy syntax, unused here because
%           we are loading pre-generated pseudoproxy data into the proxy
%           matrix instead of  generating them at runtime.
%     2) merged the hybrid RegEM and non-hybrid RegEM  cases
%           (boolean flag 'hybrid' specifies which to use)
%  02/29/12 Update (Flight to Hawaii)
%      - removed all standardization
%      - removed time reversal
%      - added version regem_2p4 (automatic choice of TTLS truncation level
%            via the MDL, AICC, or NE08 criteria)
%  11/15/12 Update 
%      - put in clean version of RegEM 2.0 as core (see regem_clean.m)
% ====================================================
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

% hybrid or not
if strmatch('hybrid', fopts)  % number of principal components retain
   hybrid = true;
   try
      f_lw = 1/options.cutoff;
   catch
      warning('No cutoff period provided for hybrid split, defaults to 20 years');
      f_lw = 0.05;
   end
else
   hybrid = false;
end

%  Principal Component Compression
if strmatch('npcs', fopts)  % number of principal components retain
   npcs = options.npcs;
else
   npcs = 0;  % Default = no compression (full field)
end
%   ====== end options processsing  =======


% Define Time Parameters
% ========================
year_i = min(t);  year_f = max(t);
%tm   = [year_f:-1:year_i]'; % IMPORTANT: Time runs in reverse.
nr   = size(proxy,2); % number of proxy records
[nt,ns] = size(field); % field dimensions

% Reverse time for both field and proxies
%target.tot = flipdim(field,1);
%prox.tot = flipdim(proxy,1);

target.tot = field;
prox.tot   = proxy;
inst.tot   = target.tot(calib,:);
ni         = size(inst.tot,1);


if strmatch('temp0', fopts)  % initial guess for temperature
   guess = 1;
   temp0 = options.temp0;
else
   guess = 0;
end

% Define instrumental target
%calib 		= ismember(tm,t(calib));
% Standardize target and proxies over the calibration period
%prox.tot = (prox.tot - repmat(mean(prox.tot(calib,:)),nt,1))./repmat(std(prox.tot(calib,:),0,1),nt,1);
%target.tot = (target.tot - repmat(mean(target.tot(calib,:)),nt,1))./repmat(std(target.tot(calib,:),0,1),nt,1);

if hybrid
   % Frequency split
   target.lw = zeros(nt,ns);
   display('Smoothing field');
   if guess
      for j=1:ns
         target.lw(:,j) = lowpass(target.tot(:,j),f_lw,1,1);
         temp0.lw(:,j) = lowpass(temp0(:,j),      f_lw,1,1);
      end
      temp0.hi = temp0.tot - temp0.lw;     
   else
      for j=1:ns
         target.lw(:,j) = lowpass(target.tot(:,j),f_lw,1,1);
      end
   end
   display('Smoothing proxies');
   for k = 1:nr
      prox.lw(:,k)   = lowpass(prox.tot(:,k),f_lw,1,1);
   end
   display('Done smoothing');
   % assign high-freq structures
   target.hi = target.tot - target.lw;
   prox.hi = prox.tot - prox.lw;
   %  Same split over instrumental period
   inst.lw  = target.lw(calib,:);
   inst.hi  = target.hi(calib,:);
end




%  2) Perform RegEM
%==================================

if hybrid
   if guess
      warning('Initial guesses not yet supported in hybrid mode')
   end
   
   if (npcs > 0)
      % Center Instrumental Data
      [inst.hic,inst.himean] = center(inst.hi);
      [inst.lwc,inst.lwmean] = center(inst.lw);
      % Singular Value Decomposition of Instrumental Data
      [Uh,Sh,Vh]      = svd(inst.hic,0); % "economy" decomposition
      instpcs.hi      = Uh(:,1:npcs);
      insteofs.hi     = Vh(:,1:npcs);
      instsingvals.hi = Sh(1:npcs,1:npcs);
      [Ul,Sl,Vl]      = svd(inst.lwc,0); % "economy" decomposition
      instpcs.lw      = Ul(:,1:npcs);
      insteofs.lw     = Vl(:,1:npcs);
      instsingvals.lw = Sl(1:npcs,1:npcs);
      % input matrices are instrumental PCs
      inst_mat.hi     = instpcs.hi;
      inst_mat.lw     = instpcs.lw;
      np = npcs;  % number of variables to reconstruct is npcs      
   else  % input matrices are the filtered field matrices
      inst_mat.hi     = inst.hi;
      inst_mat.lw     = inst.lw;
      np = ns; % number of variables to reconstruct is ns
   end
   %==========================================================
   % Millennnium Reconstruction
   %==========================================================
   
   % set Regularization parameter
   regpar = options.hi.regpar;
   switch regmethod
      case {'ttls'}  % if TTLS, need to determine truncation
         % with Mann 09 criterion
         if isempty(regpar)  % compute the regpar for hi&lw automatically
            % High-Frequency Truncation
            P = standardize([inst_mat.hi prox.hi(calib,:)]);
            [U,S,V] = svd(P,0); % do the decomposition, economy style. U's are the PCs
            highfevalues = diag(S.^2);
            highfpctvar = diag(S).^2./sum(diag(S).^2); % that gives the pct for each pc
            mbest=eigenselect(highfpctvar);
            trunc.hi=mbest;
            %
            % Low-Frequency Truncation
            P = standardize([inst_mat.lw prox.lw(calib,:)]);
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
   [nt,nr]             =   size(prox.hi);
   %	Assemble climate field/proxy data matrices
   Xh_in				     = NaN(nt,nr+np);
   Xh_in(calib,1:np)    = inst_mat.hi;  % Put in instrumental data
   Xh_in(:,np+1:np+nr) = prox.hi;  % the rest is (possibly incomplete) proxy data
   
   Xl_in			        = NaN(nt,nr+np);
   Xl_in(calib,1:np)    = inst_mat.lw;  % Put in instrumental data
   Xl_in(:,np+1:np+nr) = prox.lw;  % the rest is (possibly incomplete) proxy data
   
   % Perform Hybrid RegEM (Parallel if nprocs > 1) . Use block-based version (if options.*.block = 1)
   
   disp('Handling high-freq recon');
   [Xh, Mh, Ch, Xerr_h, Wh] = regem_2p0(Xh_in,options.hi);
   disp('Handling low-freq recon');
   [Xl, Ml, Cl, Xerr_l, Wl] = regem_2p0(Xl_in,options.lw);
   
   if npcs > 0
      %  SVD synthesis : field = U*S*V' = PCs*S*EOFs'; + add the mean taken out in lines 47,48
      recon.hi  = Xh(:,1:np)*instsingvals.hi*insteofs.hi' + repmat(inst.himean,[nt 1]);
      recon.lw  = Xl(:,1:np)*instsingvals.lw*insteofs.lw' + repmat(inst.lwmean,[nt 1]); % unscramble the scrambled egg
      recon.tot = recon.hi + recon.lw;
   else  % extract first np columns from RegEM output
      recon.hi  = Xh(:,1:np);
      recon.lw  = Xl(:,1:np);
      recon.tot = recon.hi + recon.lw;
   end
   avail        = ~isnan(Xh_in);
   
   diagn.err    = Xerr_l;
   diagn.w      = Wl;
   diagn.avail  = avail;
   
else
   if (npcs > 0)
      if guess
         warning('Initial guesses not yet supported in hybrid mode')
      end
      % Center Instrumental Data
      [instc,inst_mean] = center(inst.tot);
      % Singular Value Decomposition of Instrumental Data
      [Ui,Si,Vi]   = svd(instc,0); % "economy" decomposition, "i" stands for "instrumental data"
      instpcs      = Ui(:,1:npcs);
      insteofs     = Vi(:,1:npcs);
      instsingvals = Si(1:npcs,1:npcs);
      % input matrices are instrumental PCs
      inst_mat     = instpcs;
      np           = npcs;  % number of variables to reconstruct is npcs
   else  % input matrices are the filtered field matrices
      inst_mat     = inst.tot;
      np           = ns; % number of variables to reconstruct is ns
   end
   switch regmethod
      case {'ttls'}  % if TTLS, need to determine truncation with Mann 09 criterion
         if isempty(options.regpar)
            P = standardize([inst_mat prox.tot(calib,:)]);
            [U,S,V] = svd(P,0); % do the decomposition, economy style. U's are the PCs
            fpctvar = diag(S).^2./sum(diag(S).^2); % that gives the pct for each pc
            mbest=eigenselect(fpctvar);
            trunc=mbest;
            
            options.regpar = trunc;
         else
            trunc = options.regpar;
         end
         display(trunc)
   end
   %	Assemble climate field/proxy data matrices
   X_in    	          = NaN(nt,nr+np);
   X_in(calib,1:np)    = inst_mat;  % Put in instrumental data
   X_in(:,np+1:np+nr) = prox.tot;  % the rest is (possibly incomplete) proxy data
   
   if guess 
      X0 =   NaN(nt,nr+np);
      X0(:,1:np)    = temp0;  % Put in temperature data
      X0(:,np+1:np+nr) = prox.tot;  % the rest is (possibly incomplete) proxy data
      options.X0 = X0;
      options.C0 = cov(X0);
   end
   
   % Perform Hybrid RegEM (Parallel iif nprocs > 1). To use block-based
   % version, set options.*.block = true
   
   disp('Handling field recon');
   [X, M, C, Xerr,W] = regem_2p0(X_in,options);
   if npcs > 0
      %  SVD synthesis : field = U*S*V' = PCs*S*EOFs'; + add the mean taken out in lines 47,48
      recon.tot = X(:,1:np)*instsingvals*insteofs' + repmat(inst_mean,[nt 1]);
   else  % extract first np columns from RegEM output
      recon.tot = X(:,1:np);
   end
   avail = ~isnan(X_in);
   % Assign output data structures
   diagn.err = Xerr;
   %diagn.w = W;
   diagn.avail = avail;
end
% end loop

field_r = recon.tot;
%field_r = flipdim(recon.tot,1); % convert to standard Earthian time
end








