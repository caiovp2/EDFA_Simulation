function [overlap,R,intensity,intensityEr] = Function_Overlap(Fiber,Laser)
%% Function_Overlap Header
%c
%c                                                       ..'´`'..'´`..'´`..                                                   
%c                                                    
%c    function [overlap,R,intensity,intensityEr] 
%c                                        = Function_Overlap(Fiber,Laser);
%c
%c     Function_Overlap calculates the overlap parameter of a chosen
%c     wavelength using the parameters of the Fiber that is specified 
%c     by the user. The vector R,intensity and intensityEr permit the user
%c     to plot the Mode Field Intensity calculated in the function.
%c
%c     INPUT:
%c             Fiber             : Fiber Struct
%c             Laser             : Laser Struct
%c               
%c     OUTPUT: 
%c             overlap           : Overlap Parameter [] 
%c             R                 : Fiber Radius Vector [m]
%c             intensity         : Mode Field Intensity []
%c             intensityEr       : Mode Field contained in Erbium Region []
%c                 
%c     USES:   
%c             None by now.
%c     TODO:   
%c             Nothing yet.
%c
%c     SEE ALSO: EDFA_Input_Data
%c
%c
%c                                           by Caio M. Santos
%c                                           03/02/2016
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
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c   V                           : Normalized Frequency []
%c   u                           : u parameter    Ref. 3 []
%c   v                           : v parameter    Ref. 3 []
%c   Step                        : Radius distance Step [m]
%c   R1                          : Radius from Core to Cladd [m]
%c   R2                          : Radius from Cladd to End of Fiber [m]
%c   R                           : Complete Fiber Radius Vector [m]
%c   J0                          : Zero Order Bessel Function []
%c   J1                          : First Order Bessel Function []
%c   K0                          : Zero Order Modified Bessel Function []
%c   K1                          : First Order Modified Bessel Function []
%c   intensity1                  : Mode Field from Core to Cladd []
%c   intensity2                  : Mode Field from Cladd to End of Fiber []
%c   intensity                   : Complete Mode Field Intensity Vector [] 
%c   intensityEr                 : Mode Field contained in Erbium Region []
%c   overlap                     : Overlap Parameter []
%c   
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c   Laser.                           :Pump, Signal or ASE
%c       .FwdPower                    :Forward Power [W]
%c       .BwdPower                    :Backward Power [W]
%c       .Wavelength                  :Wavelength [m]
%c   
%c   Fiber.                           :EDF(Er-Doped Fiber) Parameters
%c       .Length                      :EDF Length [m]
%c       .Radius                      :EDF Radius [m]
%c       .CoreRadius                  :EDF Core Radius [m]
%c       .ErRadius                    :EDF Erbium Radius [m]
%c       .CoreRefractive              :Core Refractive Index []
%c       .CladdRefractive             :Cladd Refractive Index []
%c       .NA                          :Numerical Aperture []
%c       .nt                          :ion Density [m^-3]
%c       .tau                         :Metastable Lifetime [s]
%c       .Aeff                        :Effective Dopant Area [m^2]
%c       .zeta                        :Saturation Parameter [(ms)^-1]
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                           IMPORTANT NOTES                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c  1) The u and v parameters are calculated based on an approximation 
%c     mentioned in reference [3]. We must check if the normalized 
%c     frequency (V) value is in between 1 and 3 (1 <= V <= 3) to use the 
%c     approximation. If this condition is not met, there may be errors 
%c     in the calculation of the overlap parameter in this script.
%c     The warning is in the "Normalized Frequency Check" section of the
%c     code.
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Checks the number of input variables

minargs  = 1; 
maxargs  = 10;
narginchk(minargs, maxargs) % Checks if the number of inputs is between
                            % minargs and maxargs

%% Code

V = 2*pi*Fiber.CoreRadius*Fiber.NA/Laser.Wavelength; % Normalized Frequency
v = 1.1428*V - 0.9960;                               % [Ref. 3]
u = sqrt(V^2 - v^2);                                 % [Ref. 3]

%Radius Vector(R1 from core to cladd. R2 from cladd to end of fiber radius)
Step = Fiber.CoreRadius/1e4;                  % Vector Size
R1 = 0               :Step:Fiber.CoreRadius;  % Radius Vector (core-cladd)
R2 = Fiber.CoreRadius:Step:Fiber.Radius;      % Radius Vector (cladd-end)
R = [R1 R2];                                  % Radius Vector (core-end)

% Bessel Equations
J0 = besselj(0,u*R1/Fiber.CoreRadius);%Zero Order Bessel Function
J1 = besselj(1,u);                    %First Order Bessel Function
K0 = besselk(0,v*R2/Fiber.CoreRadius);%Zero Order Modified Bessel Function
K1 = besselk(1,v);                    %First Order Modified Bessel Function

% Intensity distribution of the fundamental mode (Ref. 3)

% intensity1 is for R<Fiber.CoreRadius
intensity1 = v*J0/(V*J1*Fiber.CoreRadius);  % Equation at [Ref. 3]
intensity1 = intensity1.^2;                 % Equation at [Ref. 3]
intensity1 = intensity1/pi;                 % Equation at [Ref. 3]

% intensity2 is for R>Fiber.CoreRadius
intensity2 = u*K0/(V*K1*Fiber.CoreRadius);  % Equation at [Ref. 3]
intensity2 = intensity2.^2;                 % Equation at [Ref. 3]
intensity2 = intensity2/pi;                 % Equation at [Ref. 3]

% Full intensity vector
intensity = [intensity1 intensity2];        % Full vector

% Get instensity curve where R < Fiber.beff
% Intensity from R>ErRadius set to zero
intensityEr = intensity.*(R <= Fiber.beff);  

% Calculates Overlap integral through trapz function
% The integral applies if we consider that the actual erbium ion 
% distribution is constant from r=0 to r=Fiber.beff. [Ref. 2 - Pg 144]
overlap = trapz(R,R.*intensityEr);
overlap = 2*pi*overlap;


%% Normalized Frequency Check (See "Important Notes" section in header)

% The u and v parameters are calculated based on an approximation mentioned
% in reference [3]. We must check if the normalized frequency value is in
% between 1 and 3 (1 <= V <= 3) to use the approximation. If this condition
% is not met, there may be errors in the calculation of the overlap
% parameter in this script.

% Checks the Normalized Frequency Condition
if ( V <= 1 || V >= 3 )
fprintf(['Warning: The value of the Normalized Frequency does not meet',...
              ' the necessary\n\t\t conditions for the approximations ',...
              'used. Wavelength: ',int2str(Laser.Wavelength*1e9),' nm.',...
              '\n\t\t Please see "Important Notes" section in the',...
              ' header of Function_Overlap.\n']);
end

