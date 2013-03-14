function [Xf,vra] = wavelet_filter(X,dt,Tmin,Tmax,padq,dj,s0,J1);
%WAVELET  FILTER 
%  Uses selective scale average of the wavelet transform
%   to perform timeseries filtering  (typically bandpass) 
%   It has the great advantage of removing trends provided Tmax 
%   is not too large, something which I could not achieve with fir1.m
%
%  [Xf,vra] = wavelet_filter(X,dt,Tmin,Tmax,padq,dj,s0,J1);
%  
%
%   Computes the bandpassed version of the vector X (length N),
%   
%
%   Note : Only the Morlet wavelet (k0=6) is used.
%   For other mother wavelets, the scaling must be adjusted
%   
% INPUTS: :
%-----------
%
%    X = the time series of length N.
%    dt = amount of time between each Y value, i.e. the sampling time.
%    Tmin = minimum period retained
%    Tmax = maximum period retained
%  
%   Example : for a monthly timeseries with Tmin abd Tmax expressed
%                   in years, use dt=1/12 
%
% OUTPUTS: :
% -------------
%    Xf, the filtered timeseries, of length N
%    vra, the ratio var(Xf) /var(X)
%
% OPTIONAL INPUTS:
% 
% *** Note *** setting any of the following to -1 will cause the default
%               value to be used.
%
%    PAD = if set to 1 (default ), pad time series with enough zeroes to get
%         N up to the next higher power of 2. This prevents wraparound
%         from the end of the time series to the beginning, and also
%         speeds up the FFT's used to do the wavelet transform.
%         This will not eliminate all edge effects.
%
%    DJ = the spacing between discrete scales. Default is 0.25.
%         A smaller # will give better scale resolution, but be slower to compute.
%
%    S0 = the smallest scale of the wavelet.  Default is 2*DT.
%
%    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
%
%      Example of use :
%    dt=1/12.0 ;   % define monthly resolution year unit
%   dj = 1/6;    % this will do 6 sub-octaves per octave
%   s0 = 2*dt;    % this says start at a scale of 2 months
%    j1 = 7/dj;    % this says do 7 powers-of-two with dj su
%   Tmin=7; % in years
%   Tmax=20; % in years
%   [Xf,vra] = wavelet_filter(X,dt,Tmin,Tmax,1,dj,s0,j1);
%
%
%%  Acknowledgement :  
% ``Wavelet software was provided by C. Torrence and G. Compo,
%     and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%   (this function uses their function wavelet.m, which uses wave_bases.m)
%
%   History : created 08/03/2006, Julien Emile-Geay, LDEO
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nargin < 8), J1 = -1;, end
if (nargin < 7), s0 = -1;, end
if (nargin < 6), dj = -1;, end
if (nargin < 5), padq = 1;, end
if (nargin < 4)
   error('Must input a vector X and sampling time DT, Tmin and Tmax')
end

if (size(X,2) ~=1), X = X';end
    
n1 = length(X);

if (s0 == -1), s0=2*dt;, end
if (dj == -1), dj = 1./4.;, end
if (J1 == -1), J1=fix((log(n1*dt/s0)/log(2))/dj);, end


%  compute the wavelet transform with Morlet wavelet

[wave,period,scale,coi] = wavelet(X,dt,padq,dj,s0,J1,'Morlet');
Cdelta = 0.776; % this is for the MORLET wavelet

% Scale-average between periods  Tmin Tmax 
avg = find((period >= Tmin) & (period <= Tmax));
scale_rt = sqrt((scale'))*(ones(1,n1));  % expand scale --> (J+1)x(N) array
Xs=sum(real(wave(avg,:))./scale_rt(avg,:),1); % [Eqn(29) in Torrence & Compo 1998]
%  finish
Xf = dj*sqrt(dt)/(Cdelta*pi^(-1/4))*Xs; % only valid for Morlet : change if other wavelet is used
% compute variance ratio
vra=std(Xf)^2/std(X)^2;

if size(X,1) == length(X);
	Xf = Xf';
end	
return


%WAVELET  1D Wavelet transform with optional singificance testing
%
%   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
%
%   Computes the wavelet transform of the vector Y (length N),
%   with sampling rate DT.
%
%   By default, the Morlet wavelet (k0=6) is used.
%   The wavelet basis is normalized to have total energy=1 at all scales.
%
%
% INPUTS:
%
%    Y = the time series of length N.
%    DT = amount of time between each Y value, i.e. the sampling time.
%
% OUTPUTS:
%
%    WAVE is the WAVELET transform of Y. This is a complex array
%    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
%    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
%    The WAVELET power spectrum is ABS(WAVE)^2.
%    Its units are sigma^2 (the time series variance).
%
%
% OPTIONAL INPUTS:
% 
% *** Note *** setting any of the following to -1 will cause the default
%               value to be used.
%
%    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
%         N up to the next higher power of 2. This prevents wraparound
%         from the end of the time series to the beginning, and also
%         speeds up the FFT's used to do the wavelet transform.
%         This will not eliminate all edge effects (see COI below).
%
%    DJ = the spacing between discrete scales. Default is 1/6..
%         A smaller # will give better scale resolution, but be slower to plot.
%
%    S0 = the smallest scale of the wavelet.  Default is 2*DT.
%
%    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
%
%    MOTHER = the mother wavelet function.
%             The choices are 'MORLET', 'PAUL', or 'DOG'
%
%    PARAM = the mother wavelet parameter.
%            For 'MORLET' this is k0 (wavenumber), default is 6.
%            For 'PAUL' this is m (order), default is 4.
%            For 'DOG' this is m (m-th derivative), default is 2.
%
%
% OPTIONAL OUTPUTS:
%
%    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
%           to the SCALEs.
%
%    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
%            where J1+1 is the total # of scales.
%
%    COI = if specified, then return the Cone-of-Influence, which is a vector
%        of N points that contains the maximum period of useful information
%        at that particular time.
%        Periods greater than this are subject to edge effects.
%        This can be used to plot COI lines on a contour plot by doing:
%
%              contour(time,log(period),log(power))
%              plot(time,log(coi),'k')
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%
% Notice: Please acknowledge the use of this program in any publications:
%   ``Wavelet software was provided by C. Torrence and G. Compo,
%     and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Notice: Please acknowledge the use of the above software in any publications:
%    ``Wavelet software was provided by C. Torrence and G. Compo,
%      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
%            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
%
% Please send a copy of such publications to either C. Torrence or G. Compo:
%  Dr. Christopher Torrence               Dr. Gilbert P. Compo
%  Advanced Study Program                 NOAA/CIRES Climate Diagnostics Center
%  National Center for Atmos. Research    Campus Box 216
%  P.O. Box 3000                          University of Colorado at Boulder
%  Boulder CO 80307--3000, USA.           Boulder CO 80309-0216, USA.
%  E-mail: torrence@NOSPAMucar.edu        E-mail: gpc@NOSPAMcdc.noaa.gov
%     (remove NOSPAM from email addresses)
%----------------------------------------------------------------------------
function [wave,period,scale,coi] = ...
	wavelet(Y,dt,pad,dj,s0,J1,mother,param);

if (nargin < 8), param = -1;, end
if (nargin < 7), mother = -1;, end
if (nargin < 6), J1 = -1;, end
if (nargin < 5), s0 = -1;, end
if (nargin < 4), dj = -1;, end
if (nargin < 3), pad = 0;, end
if (nargin < 2)
	error('Must input a vector Y and sampling time DT')
end

n1 = length(Y);

if (s0 == -1), s0=2*dt;, end
if (dj == -1), dj = 1./6.;, end
if (J1 == -1), J1=fix((log(n1*dt/s0)/log(2))/dj);, end
if (mother == -1), mother = 'MORLET';, end

%....construct time series to analyze, pad if necessary
x(1:n1) = Y - mean(Y);
if (pad == 1)
	base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
	x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

%....construct wavenumber array used in transform [Eqn(5)]
k = [1:fix(n/2)];
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

%....compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

%....construct SCALE array & empty PERIOD & WAVE arrays
scale = s0*2.^((0:J1)*dj);
period = scale;
wave = zeros(J1+1,n);  % define the wavelet array
wave = wave + i*wave;  % make it complex

% loop through all scales and compute transform
for a1 = 1:J1+1
	[daughter,fourier_factor,coi,dofmin]=wave_bases(mother,k,scale(a1),param);	
	wave(a1,:) = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
end

period = fourier_factor*scale;
coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]
wave = wave(:,1:n1);  % get rid of padding before returning

return


%WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG
%
%  [DAUGHTER,FOURIER_FACTOR,COI,DOFMIN] = ...
%      wave_bases(MOTHER,K,SCALE,PARAM);
%
%   Computes the wavelet function as a function of Fourier frequency,
%   used for the wavelet transform in Fourier space.
%   (This program is called automatically by WAVELET)
%
% INPUTS:
%
%    MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
%    K = a vector, the Fourier frequencies at which to calculate the wavelet
%    SCALE = a number, the wavelet scale
%    PARAM = the nondimensional parameter for the wavelet function
%
% OUTPUTS:
%
%    DAUGHTER = a vector, the wavelet function
%    FOURIER_FACTOR = the ratio of Fourier period to scale
%    COI = a number, the cone-of-influence size at the scale
%    DOFMIN = a number, degrees of freedom for each point in the wavelet power
%             (either 2 for Morlet and Paul, or 1 for the DOG)
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%----------------------------------------------------------------------------
function [daughter,fourier_factor,coi,dofmin] = ...
	wave_bases(mother,k,scale,param);

mother = upper(mother);
n = length(k);

if (strcmp(mother,'MORLET'))  %-----------------------------------  Morlet
	if (param == -1), param = 6.;, end
	k0 = param;
	expnt = -(scale.*k - k0).^2/2.*(k > 0.);
	norm = sqrt(scale*k(2))*(pi^(-0.25))*sqrt(n);    % total energy=N   [Eqn(7)]
	daughter = norm*exp(expnt);
	daughter = daughter.*(k > 0.);     % Heaviside step function
	fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h]
	coi = fourier_factor/sqrt(2);                  % Cone-of-influence [Sec.3g]
	dofmin = 2;                                    % Degrees of freedom
elseif (strcmp(mother,'PAUL'))  %--------------------------------  Paul
	if (param == -1), param = 4.;, end
	m = param;
	expnt = -(scale.*k).*(k > 0.);
	norm = sqrt(scale*k(2))*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n);
	daughter = norm*((scale.*k).^m).*exp(expnt);
	daughter = daughter.*(k > 0.);     % Heaviside step function
	fourier_factor = 4*pi/(2*m+1);
	coi = fourier_factor*sqrt(2);
	dofmin = 2;
elseif (strcmp(mother,'DOG'))  %--------------------------------  DOG
	if (param == -1), param = 2.;, end
	m = param;
	expnt = -(scale.*k).^2 ./ 2.0;
	norm = sqrt(scale*k(2)/gamma(m+0.5))*sqrt(n);
	daughter = -norm*(i^m)*((scale.*k).^m).*exp(expnt);
	fourier_factor = 2*pi*sqrt(2./(2*m+1));
	coi = fourier_factor/sqrt(2);
	dofmin = 1;
else
	error('Mother must be one of MORLET,PAUL,DOG')
end

return



