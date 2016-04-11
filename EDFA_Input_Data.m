%% EDFA Simulation Input Data
%c
%c                                                       ..'´`'..'´`..'´`..                                                   
%c     File: EDFA_Input_Data.m
%c
%c     EDFA Simulation Input Data. 
%c
%c                                           by Caio M. Santos
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
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c   h                                :Planck's Constant [Js]
%c   m                                :Number of modes []
%c   c                                :Speed of light [m/s]
%c                                    
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                                     
%c                                         
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                                    
%c   Pump.                            :Pump Laser
%c       .Power                       :Power [W]
%c       .u                           :Direction of Propagation []
%c       .l                           :Excess Losses [dB/m]
%c       .sigmaA                      :Absorption Cross Section [m^2]
%c       .sigmaE                      :Emission Cross Section [m^2]
%c       .Wavelength                  :Wavelength [m]
%c  
%c   Signal.                          :Signal Laser
%c       .Power                       :Power [W]
%c       .u                           :Direction of Propagation []
%c       .l                           :Excess Losses [dB/m]
%c       .sigmaA                      :Absorption Cross Section [m^2]
%c       .sigmaE                      :Emission Cross Section [m^2]
%c       .Wavelength                  :Wavelength [m]
%c
%c   ASE.                             :ASE
%c       .Power                       :Power [W]
%c       .BW                          :Bandwidth [Hz]
%c       .u                           :Direction of Propagation []
%c       .l                           :Excess Losses [dB/m]
%c       .sigmaA                      :Absorption Cross Section [m^2]
%c       .sigmaE                      :Emission Cross Section [m^2]
%c       .Wavelength                  :Wavelength [m]
%c   
%c   Fiber.                           :EDF(Er-Doped Fiber) Parameters
%c       .Length                      :Length [m]
%c       .Radius                      :Radius [m]
%c       .CoreRadius                  :Core Radius [m]
%c       .beff                        :Dopant Radius [m]
%c       .CoreRefractive              :Core Refractive Index []
%c       .CladdRefractive             :Cladd Refractive Index []
%c       .NA                          :Numerical Aperture []
%c       .nt                          :ion Density [m^-3]
%c       .tau                         :Metastable Lifetime [s]
%c       .Aeff                        :Effective Dopant Area [m^2]
%c       .zeta                        :Saturation Parameter [(ms)^-1]
%c   
%c   CrossSection.                    :Cross Section Data
%c       .Wavelength                  :Wavelength [m]
%c       .Absorption                  :Abs. Cross Section Spectra [m^2]      
%c       .Emission                    :Emis. Cross Section Spectra [m^2]
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                           IMPORTANT NOTES                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c

%% Global Variables

h     = 6.626070040e-34;              % Planck's Constant [Js]
c     = 299792458;                    % Speed of Light [m/s]
m     = 2;                            % Number of modes []

%% Structs

%% Signal Laser
  Signal.Power = 1e-7; % Power [W]
  Signal.u = 1;         % Direction of Propagation [1=Forward  -1=Backward]
  Signal.l = 0.03;      % Excess Losses [dB/m]
  Signal.Wavelength = 1530e-9; % Wavelength [m]
  
%% Pump Laser
  Pump.Power = 40e-3; % Power [W]
  Pump.u = 1;         % Direction of Propagation [1=Forward  -1=Backward]
  Pump.l = 0.03;      % Excess Losses [dB/m]
  Pump.Wavelength = 1480e-9; % Wavelength [m]
  
%% ASE
  ASE.Power = 0;                       % ASE Power [W]
  ASE.BW = 125e9;                       % ASE Bandwidth [Hz]
  ASE.Wavelength = [1530 1550]*1e-9;            % Wavelength [m] 
  ASE.u = 1*ones(1,length(ASE.Wavelength));            % Direction of Propagation [1=Forward  -1=Backward]
  ASE.l = 0.03*ones(1,length(ASE.Wavelength));                        % Excess Losses [dB/m]
  
%% Fiber

  % Values obtained from previous simulation (RECHECK)
  Fiber.Length = 14;;              % Length [m]                -RECHECK
  Fiber.Radius = 5e-06;;             % Radius [m]                -RECHECK
  Fiber.CoreRadius = 1.4e-06;;       % Core Radius [m]           -RECHECK
  Fiber.beff = 1.05e-06;;             % Dopant Radius [m]         -RECHECK
  Fiber.CoreRefractive = 1.45;      % Core Refractive Index []  -RECHECK
  Fiber.CladdRefractive = 1.42;;     % Cladd Refractive Index [] -RECHECK
  Fiber.NA = 0.28; %sqrt(Fiber.CoreRefractive^2-Fiber.CladdRefractive^2);;                  % Numerical Aperture []     -RECHECK
  Fiber.nt = 0.7e25;;               % ion density [m^-3]        -RECHECK
  Fiber.tau = 10e-3;;                % Metastable Lifetime [s]   -RECHECK
  Fiber.Aeff = pi*Fiber.beff^2;      % Effective Dopant Area [m^2]
  Fiber.zeta = Fiber.Aeff*Fiber.nt/Fiber.tau;;%Saturation Parameter [(ms)^-1]
  
%% Cross Section

  %load('Validation_Tests/Test_Cross_Sections\Cross.txt');
  load('Fiber_Data/FiberA/FiberA_CrossSection.txt');
  Cross = FiberA_CrossSection;
  clear FiberA_CrossSection
  
  CrossSection.Wavelength = Cross(:,1)*1e-9;    % Wavelength Vector [m]
  CrossSection.Absorption = Cross(:,2); % Absorption Cross Section [m^2]
  CrossSection.Emission   = Cross(:,3); % Emission Cross Section [m^2]
  clear Cross  
