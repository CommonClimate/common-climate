function varargout=wtc_r(x,y,varargin)
%% Wavelet coherence
%
% USAGE: [Rsq,period,scale,coi,sig95]=wtc_r(x,y,[,settings])
%
% 
% Settings: Pad: pad the time series with zeros? 
% .         Dj: Octaves per scale (default: '1/12')
% .         S0: Minimum scale
% .         J1: Total number of scales
% .         Mother: Mother wavelet (default 'morlet')
% .         MaxScale: An easier way of specifying J1
% .         MakeFigure: Make a figure or simply return the output.
% .         BlackandWhite: Create black and white figures
% .         AR1: the ar1 coefficients of the series 
% .              (default='auto' using a naive ar1 estimator. See ar1nv.m)
% .         MonteCarloCount: Number of surrogate data sets in the significance calculation. (default=1000)
% .         ts_ref: reference timeseries for Monte Carlo simulation (default = x)
% .         ArrowDensity (default: [30 30])
% .         ArrowSize (default: 1)
% .         ArrowHeadSize (default: 1)
%
% Settings can also be specified using abbreviations. e.g. ms=MaxScale.
% For detailed help on some parameters type help wavelet.
%
% Example:
%    t=1:200;
%    wtc(sin(t),sin(t.*cos(t*.01)),'ms',16)
%
% Please acknowledge the use of this software in any publications:
%   
% http://noc.ac.uk/using-science/crosswavelet-wavelet-coherence
%
% Modified Feb 21 2012 by Julien Emile-Geay (USC) to allow the specification 
% of a reference timeseries (ts_ref) against which to evaluate significance.
% See wtcsignif_r.m for the calculation
% ======================================================================





% ------validate and reformat timeseries.
[x,dt]=formatts(x);
[y,dty]=formatts(y);
if dt~=dty
    error('timestep must be equal between time series')
end
t=(max(x(1,1),y(1,1)):dt:min(x(end,1),y(end,1)))'; %common time period
if length(t)<4
    error('The two time series must overlap.')
end

ts_ref = x(:,1);

n=length(t);

%----------default arguments for the wavelet transform-----------
Args=struct('Pad',1,...      % pad the time series with zeroes (recommended)
            'Dj',1/12, ...    % this will do 12 sub-octaves per octave
            'S0',2*dt,...    % this says start at a scale of 2 years
            'J1',[],...
            'Mother','Morlet', ...
            'MaxScale',[],...   %a more simple way to specify J1
            'MakeFigure',(nargout==0),...
            'MonteCarloCount',1000,...
            'BlackandWhite',0,...
            'AR1','auto',...
            'ArrowDensity',[30 30],...
            'ArrowSize',1,...
            'ArrowHeadSize',1);
Args=parseArgs(varargin,Args,{'BlackandWhite'});
if isempty(Args.J1)
    if isempty(Args.MaxScale)
        Args.MaxScale=(n*.17)*2*dt; %auto maxscale
    end
    Args.J1=round(log2(Args.MaxScale/Args.S0)/Args.Dj);
end

ad=mean(Args.ArrowDensity);
Args.ArrowSize=Args.ArrowSize*30*.03/ad;
Args.ArrowHeadSize=Args.ArrowHeadSize*Args.ArrowSize*220;


if ~strcmpi(Args.Mother,'morlet')
    warning('Smoothing operator is designed for morlet wavelet.')
end

if strcmpi(Args.AR1,'auto')
    Args.AR1=[ar1nv(x(:,2)) ar1nv(y(:,2))];
    if any(isnan(Args.AR1))
        error('Automatic AR1 estimation failed. Specify it manually (use arcov or arburg).')
    end
end

nx=size(x,1);
sigmax=std(x(:,2));

ny=size(y,1);
sigmay=std(y(:,2));



%-----------:::::::::::::--------- ANALYZE ----------::::::::::::------------

[X,period,scale,coix] = wavelet(x(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
[Y,period,scale,coiy] = wavelet(y(:,2),dt,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);

%Smooth X and Y before truncating!  (minimize coi)
sinv=1./(scale');


sX=smoothwavelet(sinv(:,ones(1,nx)).*(abs(X).^2),dt,period,Args.Dj,scale);
sY=smoothwavelet(sinv(:,ones(1,ny)).*(abs(Y).^2),dt,period,Args.Dj,scale);


% truncate X,Y to common time interval (this is first done here so that the coi is minimized)
dte=dt*.01; %to cricumvent round off errors with fractional timesteps
idx=find((x(:,1)>=(t(1)-dte))&(x(:,1)<=(t(end)+dte)));
X=X(:,idx);
sX=sX(:,idx);
coix=coix(idx);

idx=find((y(:,1)>=(t(1))-dte)&(y(:,1)<=(t(end)+dte)));
Y=Y(:,idx);
sY=sY(:,idx);
coiy=coiy(idx);

coi=min(coix,coiy);

% -------- Cross wavelet -------
Wxy=X.*conj(Y);

% ----------------------- Wavelet coherence ---------------------------------
sWxy=smoothwavelet(sinv(:,ones(1,n)).*Wxy,dt,period,Args.Dj,scale);
Rsq=abs(sWxy).^2./(sX.*sY);

if (nargout>0)|(Args.MakeFigure)
    wtcsig=wtcsignif_r(Args.MonteCarloCount,Args.AR1,ts_ref,dt,length(t)*2,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother,.6);
    wtcsig=(wtcsig(:,2))*(ones(1,n));
    wtcsig=Rsq./wtcsig;
end

if Args.MakeFigure
    

    Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
    
    if Args.BlackandWhite
        levels = [0 0.5 0.7 0.8 0.9 1];
        [cout,H]=safecontourf(t,log2(period),Rsq,levels);

        colorbarf(cout,H)
        cmap=[0 1;.5 .9;.8 .8;.9 .6;1 .5];
        cmap=interp1(cmap(:,1),cmap(:,2),(0:.1:1)');
        cmap=cmap(:,[1 1 1]);
        colormap(cmap)
        set(gca,'YLim',log2([min(period),max(period)]), ...
            'YDir','reverse', 'layer','top', ...
            'YTick',log2(Yticks(:)), ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top')
        ylabel('Period')
        hold on

        %phase plot
        aWxy=angle(Wxy);
        lowRidx=find(Rsq<.5); %remove phase indication where Rsq is low
        aaa=aWxy;
        aaa(lowRidx)=NaN;
        %[xx,yy]=meshgrid(t(1:5:end),log2(period));

        phs_dt=round(length(t)/30); tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
        phs_dp=round(length(period)/30); pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
        phaseplot(t(tidx),log2(period(pidx)),aaa(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);

        if ~all(isnan(wtcsig))
            [c,h] = contour(t,log2(period),wtcsig,[1 1],'k');
            set(h,'linewidth',2)
        end
        %suptitle([sTitle ' coherence']);
        plot(t,log2(coi),'k','linewidth',3)
        hold off
    else
        H=imagesc(t,log2(period),Rsq);
        %[c,H]=safecontourf(t,log2(period),Rsq,[0:.05:1]);
        %set(H,'linestyle','none')
        
        set(gca,'clim',[0 1])
        
        HCB=safecolorbar;
        
        set(gca,'YLim',log2([min(period),max(period)]), ...
            'YDir','reverse', 'layer','top', ...
            'YTick',log2(Yticks(:)), ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top')
        ylabel('Period')
        hold on

        %phase plot
        aWxy=angle(Wxy);
        lowRidx=find(Rsq<.5); %remove phase indication where Rsq is low
        aaa=aWxy;
        aaa(lowRidx)=NaN;
        %[xx,yy]=meshgrid(t(1:5:end),log2(period));

        phs_dt=round(length(t)/30); tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
        phs_dp=round(length(period)/30); pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
        phaseplot(t(tidx),log2(period(pidx)),aaa(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);

        if ~all(isnan(wtcsig))
            [c,h] = contour(t,log2(period),wtcsig,[1 1],'k');
            set(h,'linewidth',2)
        end
        %suptitle([sTitle ' coherence']);
        tt=[t([1 1])-dt*.5;t;t([end end])+dt*.5];
        hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]),'w');
        set(hcoi,'alphadatamapping','direct','facealpha',.5)
        hold off
    end
end

varargout={Rsq,period,scale,coi,wtcsig};
varargout=varargout(1:nargout);






function [cout,H]=safecontourf(varargin)
vv=sscanf(version,'%i.');
if (version('-release')<14)|(vv(1)<7)
    [cout,H]=contourf(varargin{:});
else
    [cout,H]=contourf('v6',varargin{:});
end

function hcb=safecolorbar(varargin)
vv=sscanf(version,'%i.');

if (version('-release')<14)|(vv(1)<7)
    hcb=colorbar(varargin{:});
else
    hcb=colorbar('v6',varargin{:});
end

