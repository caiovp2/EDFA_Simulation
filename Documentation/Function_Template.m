function [vout1,vout2] = Function_Template(vin1,vin2,vin3,varargin);
%c
%c                                                       ..'´`'..'´`..'´`..                                                   
%c                                                    
%c    function [vout1,vout2] = Function_Template(vin1,vin2,vin3,varargin);
%c
%c     Function_Template is a template for all ONDA functions and scripts.
%c  If you are a developer, please follow this model.
%c
%c     INPUT:
%c     vin1          : Input variable 1 [Unit] 
%c             %# Please describe the variable. Do NOT forget its unit.
%c             %# Ex. vin1 : Bias voltage [Volt]
%c     vin2          : Input variable 2 [Unit] 
%c     ...
%c     vinN          : Input variable N [Unit] 
%c
%c     OUTPUT:
%c     vout1          : Output variable 1 [Unit] 
%c     vout2          : Output variable 2 [Unit] 
%c     ...
%c     voutM          : Output variable M [Unit] 
%c
%c     USES:
%c             %# Please name the ONDA functions invoked by this function.
%c     TODO: 
%c
%c     SEE ALSO: Several interesting ONDA functions
%c
%c
%c                                           by Author_Name_Surname
%c                                           Date
%c                                           email: 
%c
%c     References:
%c       [1] Segatto, M. E. V., "How to Write Nice Computer Programs: The 
%c           LabTel Example,", LabTel Press, 2015.
%c       [2] ...
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

%
% # Please, NEVER forget to check the number of input variables
% Ex.:
minargs  = 1; 
maxargs  = 10;
narginchk(minargs, maxargs) % Checks if the number of inputs is between
                            % minargs ans maxargs
%                            
% NOW you can start your CODE!
%
