%% Generate LaTeX
%c
%c                                                       ..'´`'..'´`..'´`..                                                   
%c                                                    
%c     Generates a LaTeX documentation.
%c
%c     INPUT: Document.m
%c
%c     OUTPUT: Document.tex
%c
%c     SEE ALSO: Generate_LaTeX_Documentation.m
%c
%c
%c                                           by Caio M. Santos
%c                                           20/01/2016
%c                                           caiovp2@gmail.com
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
 
%% Code

options.AuthorName = 'Caio M. Santos'; % Name of author
options.TeXheader  = 1;                % 0 means no TeX header, 1 otherwise
options.Title      = 'EDFA Simulation';% Title of document 

 nameIn  = 'EDFA.m';
 nameOut = 'EDFA.tex';

Generate_LaTeX_Documentation(nameIn,nameOut,options);

%% End

clear all