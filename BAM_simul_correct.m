function [Xc tc tmc] = BAM_simul_correct(X,t,model)
% Generate an ensemble of possible age corrected data:
% See www.clim-past-discuss.net/9/6077/2013/ for a detailed description of
% the model.
% The time series in X are automatically flipped to range from most recent to oldest measurements when the intput t is given in increasing order.
% [Xc tc tmc] = BAM_simul_correct(X,t) will generate an ensemble of 1000 age models randomly following
% a Poisson process with rate parameter theta=0.05 used to correct data X
%
% [Xc tc tmc] = BAM_simul_correct(X,t,model) will correct data X  with the model specified in
% the model structure
%
% INPUT
% X: data (vector or matrix n*p)
% t: chronology for data X (n*1)
% model.ns: number of samples
% model.name: 'poisson' or 'bernoulli'
% model.param: probability of growth band being perturbed (default: prob of missing band = prob of doubly-counted band=0.05)
%      if model.param is a single argument, then the perturbations are symmetric (prob of missing band = prob of doubly-counted band)
%      if model.param = [a1 a2] and a1 neq a2 the model is asymmetric. 
%                       a1 = prob(missing layer)
%                       a2 = prob(layer counted multiple times)
% model.resize: do not resize: 0 (default), resize to shortest sample: -1, resize to longest sample: 1
% model.tm: if a time model is provided, the code returns the corresponding corrected data

% OUTPUT
% Xc: realizations of age-perturbed data matrix of size tn*p*ns (could be 2 or 3d)
% tc: new chronology tn*1
% tmc: generated age model ensemble (tn*p*ns matrix) to map data X to the corrected realizations Xc (2=insert nan, 0=remove double band) (2 or 3d)
% where tn is the chronology length = n (default), shortest sample or longest sample
% depending on the chosen resizing option.
%
% Author: Maud Comboul, 2013. 
% Reference: http://www.clim-past-discuss.net/9/6077/2013/cpd-9-6077-2013.html
% =================================================================================

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
% Generate an ensemble of time correction models
if ~isfield(model,'tm') || isempty(model.tm)
    tmc = ones(n,p,ns);
    
    if strcmpi(model.name,'poisson')
        for nn=1:ns
            % apply poisson model
            num_event_mis = poissrnd(model.param(1)*n,[p,1]); % Poisson model for missing layers
            num_event_dbl = poissrnd(model.param(2)*n,[p,1]); % Poisson model for layers counted multiple times

            for ii=1:p
                jumps = randperm(n-1,num_event_mis(ii))+1; % place events uniformly on {2,...,n}
                tmc(jumps,ii,nn) = tmc(jumps,ii,nn)+1;   % add 1 at jump locations
                
                jumps = randperm(n-1,num_event_dbl(ii))+1; % place events uniformly on {2,...,n}
                tmc(jumps,ii,nn) = tmc(jumps,ii,nn)-1; % remove 1 at jump location
            end
            
        end
        
        
    elseif strcmpi(model.name,'bernoulli')
        for nn=1:ns
            % apply bernoulli model
            tmc(:,:,nn) = tmc(:,:,nn)+binornd(1,model.param(1),[n,p]); % Binomial model
            tmc(:,:,nn) = tmc(:,:,nn)-binornd(1,model.param(2),[n,p]); % Binomial model
        end
        
    else
        display('Unknown age model ; acceptable inputs are ''poisson'' and '' bernoulli''');
        
    end
else
    tmc = model.tm;
    ns = size(model.tm,3);
    if isflipped ==1
        for nn=1:ns
            tmc(:,:,nn) = flipud(tmc(:,:,nn));
        end
    end
end




% generate age corrected data Xc
% expand length of Xc and tXc if resizing is required
if model.resize
    t_ext = ceil(2*model.param(2)*n);
    tn = n + t_ext;
    X = [X; nan(t_ext,p)];

    dt = t(2)-t(1);
    time_ext = t(end)+dt:dt:t(end)+t_ext*dt';
    tc = [t ; time_ext'];
else
    tn=n;
    tc=t;
end

Xc = nan(tn,p,ns);
Tmax = 0;
Tmin = n;
for nn=1:ns
    for ii=1:p
        xcount=1;
        Xcount=1;
        tt=1;
        while tt<n+1
            
            if tmc(tt,ii,nn)==0     % remove double band
                Xcount=min(tn,Xcount+1);
                Xc(xcount,ii,nn)=X(min(tn,Xcount),ii);
                tt=tt+1;
            elseif tmc(tt,ii,nn)==2   % insert NaN               
                Xc(xcount,ii,nn)=nan;
                xcount=min(tn,xcount+1);
                Xc(xcount,ii,nn)=X(Xcount,ii);
            else
                Xc(xcount,ii,nn)=X(Xcount,ii);
            end
            xcount=min(tn,xcount+1);
            Xcount=min(tn,Xcount+1);
            tt=tt+1;
        end
               
        k=find(~isnan(Xc(:,ii,nn)),1,'last');
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
    Xc = Xc(1:Tmin,:,:);
    tc=tc(1:Tmin);
end
% expand output size to longest (non-nan) sequence
if model.resize==1
    Xc = Xc(1:Tmax,:,:);
    tc=tc(1:Tmax);
end

if size(X,2)==1
    Xc = reshape(Xc,n,model.ns);
    tmc = reshape(tmc,n,model.ns);
end

if isflipped == 1
    if size(X,2)==1
        Xc = flipud(Xc);
    else
        for nn=1:ns
            Xc(:,:,nn) = flipud(Xc(:,:,nn));
            tmc(:,:,nn) = flipud(tmc(:,:,nn));
        end
    end
    tc = flipud(tc);
end







