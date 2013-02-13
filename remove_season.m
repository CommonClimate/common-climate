function [anom,clim] = remove_season(array,tnum,varargin);
% function [anom,clim] = remove_season(array,tnum,[dimtime,lenseas,season_type,period]);
% 
% remove seasonal cycle from array
%
% input  -	array : data array
%        -  tnum : time axis in DATENUM FORMAT (see datenum.m)    
%		 -  dimtime	 : dimension of time [default = 1]	
%		 -  lenseas	 : length of seasonal cycle  [default = 12]
%        -  season_type : (1) 'tropical' (Apr-Mar) OR (2) 'calendar' (Jan-Dec)
%                 [default = 2]              
%        -  period : 2-element vector specifying the reference period
%                             default =  [min(year), max(year)]
%
% output -	anom : anomaly array
%		 -  clim : average seasonal cycle
% 
% History:
% =========
%  Startged by Julien Emile-Geay, USC, Oct 2011
%  Modified by N. Mirnateghi, Jan 2012 
%
%  based on remseason.m by G.Krahmann, LODYC Paris, Mar 1998  
% 
%   modified to include non-January-to-December endpoints & specific
%   reference period for the computation of the climatology. 
%   NB: assumes time runs from old to recent, unlike most paleoclimate records
%
%  Further modification: using global month variable in order to make sure
%  starting at the correct monthly cycle
%
% updated 9 Jan 2011 - nasim@earth.usc.edu
% ===================================================================

% time axis
tvec  = datevec(tnum(:));
year  = tvec(:,1);
month = tvec(:,2);

% check input
if nargin == 6
    dimtime = varargin{1};
    lenseas = varargin{2};
    season_type = varargin{3};
    period  = varargin{4};

elseif nargin == 5
    dimtime = varargin{1};
    lenseas = varargin{2};
    season_type = varargin{3};
    period  = [min(year),max(year)];      
    
elseif nargin == 4
    dimtime = varargin{1};
    lenseas = varargin{2};
    season_type = 2;
    period = [min(year),max(year)];
    
elseif nargin == 3
    dimtime = varargin{1};
    lenseas = 12;
    season_type = 2;
    period = [min(year),max(year)];
    
elseif nargin == 2
    dimtime = 1;
    lenseas = 12;
    season_type = 2;
    period = [min(year),max(year)];
else
    error('Thou shalt specify some arguments!')
end

% prepare output arrays
sd = size(array);
sn = sd;
sn(dimtime) = lenseas;
anom = repmat(0,sd);
clim = repmat(0,sn);

% shift dimensions to get the wanted one to the left
array = shiftdim(array,dimtime-1);
clim = shiftdim(clim,dimtime-1);
anom = shiftdim(anom,dimtime-1);
sds = shiftlr(sd,dimtime-1);
sns = shiftlr(sn,dimtime-1);

% calculate climatology

mRes=ceil(12/lenseas);

mgrid(1,:)=month(1);
for c = 2:lenseas
    if mod(mgrid(c-1,:)+mRes,12)==0
        mgrid(c,:)=12;
    else
        mgrid(c,:)=mod(mgrid(c-1,:)+mRes,12);
    end
end

mgrid=sort(mgrid);  % for Jan-Dec years 
 
if season_type == 1 % for May-Apr years (Tropical Years)
    while mgrid(1)<5
        mgrid=circshift(mgrid,-1);
    end
end

% fixing the start and end of the years to have a proper and complete
% seasonal cycle

month_c = month; % month_c has the complete years (Jan - Dec or Apr-Mar)
year_c  = year;  
i=1;    
while month(i) ~= mgrid(1)  
    month_c(i)=[];
    year_c(i)=[];
    i=i+1;
end

i=length(month_c);
while month_c(i)~=mgrid(length(mgrid))
    month_c(i)=[]; year_c(i)=[];
    i=i-1;
end
      
for imon = 1:lenseas 
  mon   = find(month_c == mgrid(imon));  % indices matching that particular month
  nm    = length(mon);
  mon_p = find(month_c == mgrid(imon) & ismember(year_c,[period(1):period(2)]));
  clim(imon,:) = nmean(array(mon_p,:),1);
  anom(mon,:) = array(mon,:) - repmat(clim(imon,:),[nm 1]);
end

% calculate anomalies
if isint(sd(dimtime)/lenseas)
  anom = array - repmat(clim,[sd(dimtime)/lenseas,ones(1,length(sd)-1)]);
else
  ind = [1:floor(sd(dimtime)/lenseas)*lenseas];
  res = [1:sd(dimtime)-floor(sd(dimtime)/lenseas)*lenseas];
  anom(ind,:) = array(ind,:) - ...
	repmat(clim(1:lenseas,:),...
	[floor(sd(dimtime)/lenseas),ones(1,length(sd)-1)]);
  anom(ind(length(ind))+res,:) = array(ind(length(ind))+res,:) - clim(res,:);
end

% reshape and reshift results
anom = reshape(anom,sds);
clim = reshape(clim,sns);
anom = shiftdim(anom,ndims(anom)-dimtime+1);
clim = shiftdim(clim,ndims(anom)-dimtime+1);
return
