function [Xp tp tmc] = BAM_simul_perturb(X,t,model)
% Generate an ensemble of age perturbed data.See www.clim-past-discuss.net/9/6077/2013/ for a detailed description of
% the model.
% The time series in X are automatically flipped to range from most recent to oldest measurements when the intput t is given in increasing order. 
%
% [Xp tp tmc] = BAM_simul_perturb(X,t) will generate an ensemble of 1000 age models randomly following
% a Poisson process with rate parameter theta=0.05 used to perturb data X
%
% [Xp tp tmc] = BAM_simul_perturb(X,t,model) will perturb data X  with the model specified in
% the model structure
%
% INPUT
% X: data (vector or matrix n*p)
% t: chronology for data X (n*1)
% model.ns: number of samples
% model.name: 'poisson' or 'bernoulli'
% model.param: probability of growth band being perturbed (default: prob of missing band = prob of doubly-counted band = 0.05)
%      if model.param is a single argument, then the perturbations are symmetric (prob of missing band = prob of doubly-counted band)
%      if model.param = [a1 a2] and a1 neq a2 the model is asymmetric
%                       a1 = prob(missing layer)
%                       a2 = prob(layer counted multiple times)
% model.resize: do not resize: 0 (default), resize to shortest sample: -1, resize to longest sample: 1
% model.tm: if a time model is provided, the code returns the corresponding perturbed data

% OUTPUT
% Xp: realizations of age-perturbed data matrix of size tn*p*ns (could be 2 or 3d)
% tp: new chronology tn*1
% tmc: corresponding ensemble of time-correction matrices (tn*p*ns) to map realizations in Xp back to the original data X (2=insert nan, 0=remove double band) (2 or 3d)
% where tn is the chronology length = n (default), shortest sample or longest sample
% depending on the chosen resizing option.

if size(X,1) < size(X,2)
    X = X';
end

t=t(:);
if mean(diff(t))>0
    isflipped = 1;
    X = flipud(X);
    t = flipud(t);
else
    isflipped = 0;
end
[n p]=size(X);

if nargin < 3
    model.name='poisson';
    model.param=[0.05 0.05];
    model.ns=1000;
end

if ~isfield(model,'ns')
    model.ns=1000;
end

if ~isfield(model,'name')
        model.name='poisson';
end

if ~isfield(model,'param')
        model.param=[0.05 0.05];
end

if ~isfield(model,'resize')
    model.resize=0;
end

if length(model.param)==1
    model.param(2) = model.param(1);
end

ns=model.ns;

% Generate an ensemble of time perturbation models
if ~isfield(model,'tm') || isempty(model.tm)
    model.tm = ones(n,p,ns);
    if strcmpi(model.name,'poisson')
        % apply poisson model
        for nn=1:ns
            
            num_event_mis = poissrnd(model.param(1)*n,[p,1]); % poisson model for missing bands 
            num_event_dbl = poissrnd(model.param(2)*n,[p,1]); % poisson model for layers counted multiple times
            
            for ii=1:p
                jumps = randperm(n-1,num_event_mis(ii))+1; % place events uniformly on {2,...,n}
                model.tm(jumps,ii,nn) = model.tm(jumps,ii,nn)-1;   % remove 1 at jump locations
                jumps = randperm(n-1,num_event_dbl(ii))+1;
                model.tm(jumps,ii,nn) = model.tm(jumps,ii,nn)+1;   % add 1 at jump locations
            end
            
        end
        
    elseif strcmpi(model.name,'bernoulli')
        for nn=1:ns
            
            % apply bernoulli model
            model.tm = model.tm-binornd(1,model.param(1),[n,p]); % Binomial model for missing bands 
            model.tm = model.tm+binornd(1,model.param(2),[n,p]); % Binomial model for layers counted multiple times
            
        end
        
    else
        display('Unknown age model ; acceptable inputs are ''poisson'' and '' bernoulli''');
        
    end
end



% generate age perturbed data Xp
% expand length of Xp and tXp if resizing is required
if model.resize
    t_ext = ceil(2*model.param(2)*n);
    tn = n + t_ext;
    X = [X; nan(t_ext,p)];

    dt = t(2)-t(1);
    time_ext = t(end)+dt:dt:t(end)+t_ext*dt';
    tp = [t ; time_ext'];
else
    tn=n;
    tp=t;
end

Xp = nan(tn,p,ns);
Tmax = 0;
Tmin = n;

tmc=ones(tn,p,ns);
for nn=1:ns
    
    for ii=1:p
        xcount=1;
        Xcount=1;
        tt=1;
        while tt < n+1
            
            if model.tm(tt,ii,nn)==0     % remove band
                Xcount=min(Xcount+1,tn);
                tmc(xcount,ii,nn)=tmc(xcount,ii,nn)+1;
            elseif   model.tm(tt,ii,nn)==2   % insert double band
                Xp(xcount,ii,nn)=X(Xcount,ii);
                tmc(xcount,ii,nn)=tmc(xcount,ii,nn)-1;
                xcount=min(tn,xcount+1);
            end
            
            Xp(xcount,ii,nn)=X(Xcount,ii);
            xcount=min(tn,xcount+1);
            Xcount=min(tn,Xcount+1);
            tt=tt+1;
        end
        k=find(~isnan(Xp(:,ii,nn)),1,'last');
        if k > Tmax
            Tmax=k;
        end
        if k < Tmin
            Tmin=k;
        end
    end
end
% crop output size to shortest sequence
if model.resize==-1   
    n=Tmin;
end
% expand output size to longest (non-nan) sequence
if model.resize==1
    n=Tmax;
end
Xp = Xp(1:n,:,:);
tp=tp(1:n);
tmc = tmc(1:n,:,:);
    
if size(X,2)==1
    Xp = reshape(Xp,n,model.ns);
    tmc = reshape(tmc,n,model.ns);
end

if isflipped == 1
    if size(X,2)==1
        Xp = flipud(Xp);
    else
        for nn=1:ns
            Xp(:,:,nn) = flipud(Xp(:,:,nn));
            tmc(:,:,nn) = flipud(tmc(:,:,nn));
        end
    end
    tp = flipud(tp);
end

