%% EDFA Simulation
%c
%c                                                       ..'´`'..'´`..'´`..                                                   
%c     File: EDFA.m
%c
%c     Simulation of Erbium-Doped Fiber Amplifiers.
%c
%c     INPUT:
%c
%c     OUTPUT: 
%c
%c     SEE ALSO: EDFA_Input_Data
%c
%c
%c                                           Author: Caio M. Santos
%c                                           20/01/2016
%c                                           caiovp2@gmail.com
%c
%c     References:
%c       [1] Freitas, M. , "Amplificadores Óticos a Fibra sob um Ambiente 
%c           Dinâmico", LabTel Press, 2006.
%c       [2] P. C. Becker, N.A. Olsson, and J. R. Simpson. "Erbium-Doped
%c           Fiber Amplifiers: Fundamentals and Technology". Optics and
%c           Photonics, 1999.
%c       [3] C. Randy Giles, and E. Desurvire. "Modeling Erbium-Doped 
%c           Fiber Amplifiers". IEEE Journal of Lightwave Technology,
%c           Volume 9, No. 2, Feb. 1991, pp. 271-283.
%c       [4] C. R. Giles, C.A. Burrus, D.J. DiGiovanni, N.K. Dutta, and
%c           G. Raybon. "Characterization of Erbium-Doped Fibers and 
%c           Application to Modeling 980 nm and 1480 nm Pumped Amplifiers".
%c           IEEE Photonics Technology Letters, Volume 3, No. 4, Apr. 1991,
%c           pp. 363-365.
%c           
%c
%c    This file is part of ONDA, the Optical Network Design and Analysis 
%c    tool.
%c    Copyright (C) 2015  LabTel - Laboratorio de Telecomunicacoes
%c                   Federal University of Espirito Santo - Brazil
%c                                   http://www.labtel.ele.ufes.br
%c                                             segatto@ele.ufes.br
%c			 
%c    ONDA is free software; you can redistribute it and/or modify
%c    it under the terms of the GNU General Public License as published by
%c    the Free Software Foundation; either version 3 of the License, or
%c    (at your option) any later version.
%c
%c    ONDA is distributed in the hope that it will be useful,
%c    but WITHOUT ANY WARRANTY; without even the implied warranty of
%c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%c    GNU General Public License for more details.
%c
%c    You should have received a copy of the GNU General Public License
%c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
                            
%% Input Files

% Starts counter to time the execution
tic;

% Wipes out data from previous executions
clear all; % Cleans memory
clc;       % Cleans screen
close all; % Closes figures

% Calls Inputs Files
EDFA_Input_Data;

%% Calculates the Overlap Parameters

for kk=1:length(Pump.Wavelength)
[Pump.Overlap(:,kk),Pump.RC(:,kk),Pump.MI(:,kk),~] = Function_Overlap(Fiber,Pump,kk);
end

for kk=1:length(Signal.Wavelength) 
[Signal.Overlap(:,kk),Signal.RC(:,kk),Signal.MI(:,kk),~] = Function_Overlap(Fiber,Signal,kk);
end

for kk=1:length(ASE.Wavelength)
[ASE.Overlap(:,kk),~,~,~] = Function_Overlap(Fiber,ASE,kk);
end

clear kk

%% Absorption and Gain coefficients

for kk=1:length(Pump.Wavelength)
Pump.Absorption(1,kk)   = Pump.Overlap(1,kk)*Pump.sigmaA(1,kk)*Fiber.nt;
Pump.Gain(1,kk)         = Pump.Overlap(1,kk)*Pump.sigmaE(1,kk)*Fiber.nt;
end

for kk=1:length(Signal.Wavelength)
Signal.Absorption(1,kk) = Signal.Overlap(1,kk)*Signal.sigmaA(1,kk)*Fiber.nt;
Signal.Gain(1,kk)       = Signal.Overlap(1,kk)*Signal.sigmaE(1,kk)*Fiber.nt;
end

for kk=1:length(ASE.Wavelength)
ASE.Absorption(1,kk) = ASE.Overlap(1,kk)*ASE.sigmaA(1,kk)*Fiber.nt;
ASE.Gain(1,kk)       = ASE.Overlap(1,kk)*ASE.sigmaE(1,kk)*Fiber.nt;
end

clear kk

%% Analytical Solution to the two-level system

% Qin = [Signal.Power*Signal.Wavelength Pump.Power*Pump.Wavelength]./(h*c);
% 
% implicit = @(x,y) implicit_func(x,Fiber,Signal,Pump,h,c);
% xo = [0 sum(Qin)];
% Qout = fzero(implicit,xo);
% 
% A = [Signal.Absorption Pump.Absorption];
% G = [Signal.Gain Pump.Gain];
% 
% Ps = Qin(1,1)*exp( ( A(1,1)+G(1,1) )*( Qin(1,1)+Qin(1,2)-Qout )/Fiber.zeta -A(1,1)*Fiber.Length )*h*c/Signal.Wavelength; 
% Pp = Qin(1,2)*exp( ( A(1,2)+G(1,2) )*( Qin(1,1)+Qin(1,2)-Qout )/Fiber.zeta -A(1,2)*Fiber.Length )*h*c/Pump.Wavelength;
% 
% clear Qin Qout xo A G

%% Solving Differential Equations
aux = 1;
for ss=1:length(Signal.Wavelength);    
for pp=1:length(Pump.Wavelength);

    for kk=1:4
        Pump.Power = kk*10e-3;
    odes = @(x,y) diff_func(x,y,Fiber,Signal,Pump,ASE,h,m,c,ss,pp);
    bcs = @(xa,ya) bcs_func(xa,ya,Signal,Pump,ASE);

    %options = bvpset('RelTol',1e-9,'AbsTol',1e-9,'Nmax',1e3);
    options = bvpset('RelTol',1e-5,'AbsTol',1e-9);

    %solinit=bvpinit(linspace(0,Fiber.Length,50),[ 10e-3 1e-4 0 0]);
    solinit=bvpinit(linspace(0,Fiber.Length,50),[ Pump.Power/2 Signal.Power 1e-3 1e-3]);

    sol = bvp4c(odes,bcs,solinit,options);
    gain(kk,aux) = 10*log10(sol.y(1,length(sol.y))/sol.y(1,1));
    
    
         N2 =((h*Fiber.zeta)^-1)*(sol.y(:,1)*Signal.Absorption(1,ss)/(c/Signal.Wavelength(1,ss))+...
                              sol.y(:,2)*Pump.Absorption(1,pp)/(c/Pump.Wavelength(1,pp))+...
                              sol.y(:,3)*ASE.Absorption(1,ss)/(c/ASE.Wavelength(1,ss))+...
                           sol.y(:,4)*ASE.Absorption(1,ss)/(c/ASE.Wavelength(1,ss)))*...
                            Fiber.nt./(1 +((h*Fiber.zeta)^-1)*...
                (sol.y(:,1)*(Signal.Absorption(1,ss)+Signal.Gain(1,ss))/(c/Signal.Wavelength(1,ss))+...
                sol.y(:,2)*(Pump.Absorption(1,pp)+Pump.Gain(1,pp))/(c/Pump.Wavelength(1,pp))+...
                sol.y(:,3)*(ASE.Absorption(1,ss)+ASE.Gain(1,ss))/(c/ASE.Wavelength(1,ss))+...
                sol.y(:,4)*(ASE.Absorption(1,ss)+ASE.Gain(1,ss))/(c/ASE.Wavelength(1,ss))));
                gain2(kk,aux) = exp(((Signal.Absorption(1,ss)+Signal.Gain(1,ss))*mean(N2)-Signal.Absorption(1,ss))*Fiber.Length);
    
    end
aux = aux + 1;
end
end

plot(gain)
legend('980:1530','1480:1530','980:1550','1480:1550')

%% Later ( COMENTADO ) (NOT FINISHED)

% N2 =((h*Fiber.zeta)^-1)*(Power(:,1)*Signal.Absorption/Signal.Wavelength+...
%                          Power(:,2)*Pump.Absorption/Pump.Wavelength+...
%                          Power(:,3)*ASE.Absorption/ASE.Wavelength+...
%                          Power(:,4)*ASE.Absorption/ASE.Wavelength)*...
%       Fiber.nt./(1 +((h*Fiber.zeta)^-1)*...
%      (Power(:,1)*(Signal.Absorption+Signal.Gain)/Signal.Wavelength+...
%       Power(:,2)*(Pump.Absorption+Pump.Gain)/Pump.Wavelength+...
%       Power(:,3)*(ASE.Absorption+ASE.Gain)/ASE.Wavelength+...
%       Power(:,4)*(ASE.Absorption+ASE.Gain)/ASE.Wavelength));
%   
% gain = 10*log10(Power(:,1)/Power(1,1));
% G = exp((Signal.Absorption+Signal.Gain)*mean(N2)/Fiber.nt-Signal.Absorption)*Fiber.Length;

%% End

% Display execution time
fprintf('\nTempo de processamento: %d segundos.\n',floor(toc));

return; % Ends the simulation

% Below are some figures for analysing separately. 
% Use Ctrl+Enter to execute each section and see the desired figure.

%% [Figure 1] Cross Sections

figure;
plot(CrossSection.Wavelength*1e9,CrossSection.Absorption,'k',...
     CrossSection.Wavelength*1e9,CrossSection.Emission,'--k',...
     'LineWidth',1.5);
axis([900 1650 0 10e-25]); 
legend('Absorption','Emission')
title('Absorption and Emission Cross Sections')
xlabel('Wavelength [nm]')
ylabel('Cross Section [m^2]')
set(gca,'YTick',0:2e-25:10e-25)

%% [Figure 2] Mode Intensity

figure;
plot(Pump.RC*1e6,Pump.MI,'LineWidth',2);
hold on;
plot(Signal.RC*1e6,Signal.MI,'r','LineWidth',2);
if(Fiber.beff~=Fiber.CoreRadius)
plot(Fiber.beff*1e6*ones(1,2),[0 max([Pump.MI Signal.MI])],'-.k',...
     'LineWidth',2);
plot(Fiber.CoreRadius*1e6*ones(1,2),[0 max([Pump.MI Signal.MI])],':k',...
     'LineWidth',2);
legend('Pump','Signal','Confined erbium boundary','Core boundary');
else
plot(Fiber.CoreRadius*1e6*ones(1,2),[0 max([Pump.MI Signal.MI])],'-.k',...
     'LineWidth',2);
legend('Pump','Signal','Core and erbium boundary')
end    
xlabel('Radial Coordinate [\mum]','interpreter','Tex');
ylabel('Mode intensity');
title('Signal and Pump Mode Intensities.');
axis([0 1e6*max([Pump.RC Signal.RC]) 0 1.1*max([Pump.MI Signal.MI])]);

%% [Figure 3] Power Propagation

figure;
subplot(2,1,1);
plot(Z,Power(:,2)*1e3,'k','LineWidth',1.5);
title('Pump Propagation')
axis([0 Fiber.Length 0 1.2e3*max(Power(:,2))]);
grid on;
ylabel('Power [mW]')

subplot(2,1,2);
plot(Z,Power(:,1)*1e3,'k','LineWidth',1.5);
title('Signal Propagation')
axis([0 Fiber.Length 0 1.2e3*max(Power(:,1))]);
grid on;
xlabel('Fiber Length [m]');
ylabel('Power [mW]')

%% [Figure 4] ASE Propagation

figure;
plot(Z,Power(:,3),'k',Z,Power(:,4),'--k','LineWidth',2);
title('ASE Propagation')
legend('Forward ASE','Backward ASE');
axis([0 Fiber.Length 0 1.1*max(max(Power(:,3),Power(:,4)))]);
xlabel('Fiber Length [m]');
ylabel('Power [mW]')

%% [Figure 5] Signal Gain x Fiber Length

figure;
plot(Z,gain,'k','LineWidth',2);
grid on;
title('Signal Gain');
xlabel('Fiber Length [m]');
ylabel('Signal Gain [dB]');
axis([0 Fiber.Length 0 1.2*max(gain)])

%% [Figure 6] Population x Fiber Length

figure;
plot(Z,N2/Fiber.nt,'k','LineWidth',2);
grid on;
title('Upper State Population');
xlabel('Fiber Length [m]');
ylabel('Upper State Population');
axis([0 Fiber.Length 0 1.1]);
