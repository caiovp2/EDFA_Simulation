function varargout = receiver_cohmix(ich,x)

%RECEIVER_COHMIX Complete COHerent MIXer receiver. (POST fiber+OBPF+MIX+PD+LPF).
%   VARARGOUT=RECEIVER_COHMIX(ICH,X) returns the received current of a M-PSK 
%   transmission using the following receiver:
%
%                                                
%                        
%                                                BALANCED
%                                                   PDS
%        __            
%       /  \                         ----------    __|__     ------------------
%      |    |  -------    -----     |          |--  /_\  ---|                  |
%  sig  \  /  |       |  |     |----| COHERENT |     |      |                  |
%    ---------| OBPF  |--| PBS |  | |          |            |                  |-- amplitudes
%       post  |       |  |     | -|-|   MIXER  |   __|__    |      SIGNAL      |
%       fiber  -------    ----- | | |          |--  /_\  ---|                  |
%                               | |  ----------      |      |                  |
%                               | |                         |                  |
%                               | |                         |                  |
%                               | |  ----------    __|__    |                  |
%       -------                 | | |          |--  /_\  ---|    PROCESSING    |
%      |       |                |  -| COHERENT |     |      |                  |
%      |  LO   |                |   |          |            |                  |-- phases
%      |       |--------------------|   MIXER  |   __|__    |                  |
%       -------                     |          |--  /_\  ---|                  |
%                                    ----------      |       ------------------
%
%
%   X is a structure of fields:
%
%   X.oftype = optical filter (OBPF) type (see MYFILTER)
%   X.obw = OBPF 3 dB bandwidth normalized to the bit rate. 
%   X.oord = optical filter order (for special filters, see MYFILTER)
%   X.eftype = electrical filter (LPF) type (see MYFILTER)
%   X.ebw = LPF 3-dB bandwidth normalized to the bit rate. 
%   X.eord = electrical filter order (for special filters, see MYFILTER)
%
%   Optional parameters of X:
%
%   X.lodetuning    = Local Oscillator Detuning frequency [Hz]
%   X.lophasenoise  = Local Oscillator Phase Noise Vector [rad]
%   X.lolinewidth   = Local Oscillator Linewidth (Hz/symbolrate)
%   X.lopower       = Local Oscillator Power [dBm]
%   X.dpost  = post compensating fiber cumulated dispersion [ps/nm]
%   X.slopez = post compensating fiber cumulated slope [ps/nm^2]
%   X.lambda = wavelength [nm] at which the post compensating fiber
%              has a cumulated dispersion equal to X.dpost.
%   X.b2b = 'b2b' evaluates the currents in back-to-back configuration, i.e.
%            with the transmitter connected directly to the receiver.  
%
%   VARARGOUT has a variable number of arguments, from 1 to 2. In the
%   complete case it is VARARGOUT=[IRIC,X], where IRIC is a matrix
%   containing the received currents of channel ICH, while X is updated.
%
%   The post-compensating fiber is assumed as a purely ideal-linear fiber,
%   while the photodiode is ideal (abs(.)^2).
%
%   Note: This function works over a copy of the electric field. All fields 
%   of the global variable GSTATE are left unchanged.
%
%   See also RECEIVER_OOK, DSP4COHDEC
%
%   Author: Massimiliano Salsi, 2009
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2009  Paolo Serena, <serena@tlc.unipr.it>
%			 
%    Optilux is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    Optilux is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

global CONSTANTS;  % CONSTANTS is a global structure variable.
CLIGHT = CONSTANTS.CLIGHT;      % speed of lightin vacuum [m/s]

global GSTATE   % GSTATE is a structure whose fields are defined in reset_all.m

if ~isfield(x,'oord')
    x.oord = 0;     % not using special filter
end

if ~isfield(x,'eord')
    x.eord = 0;     % not using special filter
end

%%%%%%%%%%% INITIALIZATION

Nfft = length(GSTATE.FN);
[nfr,nfc] = size(GSTATE.FIELDX);    % if nfc == 1 there is only one field 
npoints = 1:Nfft;                   % for all channels.
maxl=max(GSTATE.LAMBDA);
minl=min(GSTATE.LAMBDA);
lamc = 2*maxl*minl/(maxl+minl); %central wavelength: 1/lamc = 0.5(1/maxl+1/minl)
if nfc ~= GSTATE.NCH
    minfreq = GSTATE.FN(2)-GSTATE.FN(1);  
    deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich)); % frequency spacing    
    ndfn = round(deltafn./GSTATE.SYMBOLRATE/minfreq);  % spacing in points
    nind = nmod(npoints-ndfn,Nfft);    % left-circular shift    
    nch = 1;
    if ich == 1 % ndfnl & ndfnr are used for evaluating the energy after
        ndfnl = Nfft/2;
    else
        deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich-1)); % left channel
        ndfnl = round(deltafn./GSTATE.SYMBOLRATE/minfreq);
        ndfnl = round((ndfn-ndfnl)*0.5);    % spacing from left chan
    end
    if ich == GSTATE.NCH
        ndfnr = Nfft/2;
    else
        deltafn = CLIGHT*(1/lamc-1/GSTATE.LAMBDA(ich+1)); % left channel
        ndfnr = round(deltafn./GSTATE.SYMBOLRATE/minfreq);        
        ndfnr = round((ndfnr-ndfn)*0.5);    % spacing from right chan
    end       
else
    nind = npoints;
    nch = ich;
    ndfnl = Nfft/2;
    ndfnr = Nfft/2;
end    
  
if isfield(x,'b2b')
    if strcmp(x.b2b,'b2b')
        b2b = 1;
        if isfield(x,'dpost'), x=rmfield(x,'dpost'); end;
    else
        error('the b2b field must be ''b2b''');
    end
else
    b2b = 0;
end

if b2b
    x.sigx = GSTATE.FIELDX_TX(:,nch);
else
    x.sigx = GSTATE.FIELDX(:,nch);
end

if isfield(x,'dpost')
    b20z = -x.lambda^2/2/pi/CLIGHT*x.dpost*1e-3; % beta2 [ns^2] @ lambda
    b30z = (x.lambda/2/pi/CLIGHT)^2*(2*x.lambda*x.dpost+...
        x.lambda^2*x.slopez)*1e-3;  
                                     % beta3 [ns^3] @ lambda

    % Domega_ik: [1/ns]. "i" -> at ch. i, "0" -> at lambda
    Domega_i0 = 2*pi*CLIGHT*(1./GSTATE.LAMBDA(ich)-1/x.lambda);    
    Domega_ic = 2*pi*CLIGHT*(1./GSTATE.LAMBDA(ich)-1/lamc);  
    Domega_c0 = 2*pi*CLIGHT*(1./lamc-1/x.lambda); 
    beta1z = b20z*Domega_ic+0.5*b30z*(Domega_i0^2-Domega_c0^2);    %[ns]    
    beta2z = b20z+b30z*Domega_i0;  % beta2*z [ns^2]@ GSTATE.LAMBDA
    % dispersion of the channels
    omega = 2*pi*GSTATE.SYMBOLRATE*GSTATE.FN';     % angular frequency [rad/ns]
    betat = omega*beta1z+0.5*omega.^2*beta2z+omega.^3*b30z/6;
    x.post_delay = GSTATE.SYMBOLRATE.*beta1z;

    Hf = fastexp(-betat);
else
    Hf = ones(Nfft,1);
    x.post_delay = 0;
end
Hf = Hf .* myfilter(x.oftype,GSTATE.FN,0.5*x.obw,x.oord);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% RECEIVER SIDE

x.sigx = fft(x.sigx);
x.sigx = x.sigx(nind);

x.avgebx = sum(abs(x.sigx(1:ndfnl)).^2) + sum(abs(x.sigx(Nfft-ndfnr+1:Nfft)).^2);
x.avgebx = x.avgebx/GSTATE.POWER(ich)/Nfft^2; % normalized average energy x bit
x.sigx = x.sigx .* Hf;                        % <- optical filter

%%%%%%%%% Calculating the local oscillator signal
% We have three parameters: power, detuning (frequency relative to this channel
% central frequency), and phase noise. Anyone of these parameters can be optional
% and the default values are:
% Power == 0 dBm
% Detuning == 0 Hz
% Phase Noise == 0 rad  ( vector! )
if isfield(x,'lodetuning') && x.lodetuning
    %step=1/GSTATE.SYMBOLRATE/GSTATE.NT/1E9;
    %time = (0:step:(GSTATE.NSYMB*GSTATE.NT - 1)*step )';
    minfreq = GSTATE.SYMBOLRATE*1E9 / GSTATE.NSYMB;
    kdet = floor( x.lodetuning / minfreq );
    if ~kdet
        warning('optilux:receiver_cohmix',...
            'Detuning is neglected! Minimum frequency too high.');
    end
    LO_Detuning = 2*pi*kdet/Nfft*(1:Nfft)';
else
    LO_Detuning = 0;
end
if isfield(x,'lophasenoise')
    if length( x.lophasenoise ) ~= Nfft
        error('Incompatible vector.');
    end
    LO_PhaseNoise = x.lophasenoise;
elseif  isfield(x,'lolinewidth')
    freq_noise  = (ones(Nfft,1) * sqrt(2*pi*x.lolinewidth./GSTATE.NT)) .* randn( Nfft, 1);
    freq_noise(1) = 0;
    LO_PhaseNoise = cumsum(freq_noise,1);
    % Brownian bridge
    for nnoise=1:length(LO_PhaseNoise)
        LO_PhaseNoise(nnoise)=LO_PhaseNoise(nnoise)-(nnoise-1)/(length(LO_PhaseNoise)-1)*LO_PhaseNoise(end);
    end
else
    LO_PhaseNoise = zeros( Nfft, 1);
end
if isfield(x,'lopower')
    LO_Ecw = 10^( x.lopower / 20 );
else
    LO_Ecw = 1;
end

Elo = LO_Ecw * fastexp( LO_Detuning + LO_PhaseNoise); % <- Local Oscillator

isy = ~isempty(GSTATE.FIELDY);
if isy
    if b2b
        if  isempty(GSTATE.FIELDY_TX)
            x.sigy = zeros(nfr,1);
        else
            x.sigy = fft(GSTATE.FIELDY_TX(:,nch));
        end
    else
        x.sigy = fft(GSTATE.FIELDY(:,nch));
    end
    x.sigy = x.sigy(nind);
    x.avgeby = sum(abs(x.sigy(1:ndfnl)).^2) + sum(abs(x.sigy(Nfft-ndfnr+1:Nfft)).^2);
    x.avgeby = x.avgeby/GSTATE.POWER(ich)/Nfft^2; % normalized average energy x bit    
    x.sigy = x.sigy .* Hf;    
    if nargout == 2    
        varargout(2) = {x};
    end
    x.sigy = ifft(x.sigy);
    x.sigx = ifft(x.sigx);
else
    if nargout == 2    
        varargout(2) = {x};
    end
    x.sigx = ifft(x.sigx);
end

% Following the report there are 4 output fields from the coherent mixers:
% 
% MIXER POLAR X AND MIXER POLAR Y:
% 1 : Erx * exp( j* pi/2 ) + Elo * exp( j* pi/2 ) 
% 2 : Erx * exp( j* 0    ) + Elo * exp( j* pi   ) 
% 3 : Erx * exp( j* pi/2 ) + Elo * exp( j* pi   ) 
% 4 : Erx * exp( j* pi   ) + Elo * exp( j* pi/2 ) 
% where Erx is respectively sigx and sigy
EmixOutX = [x.sigx * j + Elo * j ...
            x.sigx     - Elo      ...
            x.sigx * j - Elo      ...
           -x.sigx     + Elo * j];
% The detection can be made using two balanced photodiodes or two normal
% photodiodes. If nothing is specified balanced detection is used.
% In the direct case the currents given by the fileds 1 and 3 are used.
% In the balanced case the four currents are used.
if isfield(x,'pdtype') && strcmp( x.pdtype,'normal')
    balanced = false;
else
    balanced = true;
end
IricX = real( EmixOutX .* conj(EmixOutX) ); % <- Photocurrents

if balanced
    IricX = [ IricX(:,1) - IricX(:,2) , IricX(:,3) - IricX(:,4) ];
else
    IricX = [ IricX(:,1) , IricX(:,3) ];
end

if isy      % Do the same things for the Y polar
    EmixOutY = [ x.sigy .* j + Elo .* j ...
                 x.sigy      - Elo      ...
                 x.sigy .* j - Elo      ...
                -x.sigy      + Elo .* j    ];
    IricY = real( EmixOutY .* conj(EmixOutY) );
    if balanced
        IricY = [ IricY(:,1) - IricY(:,2) , IricY(:,3) - IricY(:,4) ];
    else
        IricY = [ IricY(:,1) , IricY(:,3) ];
    end
end

% lowpass filter
Hf = myfilter(x.eftype,GSTATE.FN,x.ebw,x.eord) * [1 1]; 
% while the 3 dB bandwidth of the OBPF goes from -x.obw/2 -> +x.obw/2, 
% the LPF goes from 0 -> x.eftype.

IricX = real(ifft(fft(IricX) .* Hf));


if isy
    IricY = real(ifft(fft(IricY) .* Hf));
    varargout(1) = {[IricX IricY]};
else
    varargout(1) = {IricX};
end

