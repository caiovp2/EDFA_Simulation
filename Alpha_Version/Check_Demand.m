function out = Check_Demand(Demand);
%
%                                                        ..'´`'..'´`..'´`..                                                   
%                                                    
%     function out = Check_Demand(Demand);
%
%     Checks the struct Demand
%
%     INPUT:
%     Demand        : Demand to be checked
%            .id
%            .id_Node_Source
%            .id_Node_Destination
%            .type           : IMDD, Coherent, EOFDM, OOFDM or FlexGrid
%            .value          : in Gbit/s
%
%     OUTPUT:
%     out           : True if everything is ok, otherwise False 
%
%     USES:
%
%                                           by M.Segatto
%                                           16/04/2015
%                                           segatto@ele.ufes.br
%
%    This file is part of ONDA, the Optical Network Design and Analysis 
%    tool.
%    Copyright (C) 2015  LabTel - Laboratorio de Telecomunicacoes
%                   Federal University of Espirito Santo - Brazil
%                                   http://www.labtel.ele.ufes.br
%                                             segatto@ele.ufes.br
%			 
%    ONDA is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    ONDA is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


usar fieldnames(options);
