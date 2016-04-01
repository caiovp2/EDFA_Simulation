function Stats_Out = Generate_LaTeX_Documentation(FileIn,FileOut,options);
%c
%c                                                       ..'�`'..'�`..'�`.. 
%c                                                           .~..~~..~~..~.
%c               
%c                                                  .-.-.   .-.-.   .-.-.     
%c                                                 / / \ \ / / \ \ / / \ \  
%c                                                .-.   .-.-.   .-.-.   .-.

%c    function 
%c    Stats_Out = Generate_LaTeX_Documentation(FileIn,FileOut,options);
%c
%c     Simple LaTeX generator from  .m files
%c
%c     INPUT:
%c     FileIn        : Input m file
%c     FileOut       : Output TeX file
%c     options
%c       .AuthorName : Report author
%c       .TeXheader  : 1 - Insert TeX header; 0 - No
%c       .Title      : Report title. If not present, the input filename 
%c                     is used.
%c
%c     OUTPUT:
%c     Stats_Out     : Struct containing file statistics       
%c
%c     USES:
%c     nothing
%c
%c     TODO: 
%c
%c                                          by M.Segatto
%c                                           17/04/2015
%c                                           segatto@ele.ufes.br
%c
%c    This file is part of ONDA, the Optical Network Design and Analysis 
%c    tool.
%c    Copyright (C) 2015  LabTel - Laboratorio de Telecomunicacoes
%c                   Federal University of Espirito Santo - Brazil
%c                                   http://www.labtel.ele.ufes.br
%c                                             segatto@ele.ufes.br
%c			 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK AND INIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TeXheader = 0;
if nargin == 1
    if strfind(FileIn,'.m')
        FileOut = strcat(FileIn(1:end-2),'.tex');
    end
elseif nargin == 3
    if isfield(options,'TeXheader')
        TeXheader = options.TeXheader;
    end
    if isfield(options,'AuthorName')
        AuthorName = options.AuthorName;
    else
        Author = 'no author';
    end
    if isfield(options,'Title')
        Title = options.Title;
    else
        Title = FileIn;
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%% NOW WE CAN START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
cont     = 0;
fidIn    = fopen(FileIn);
fidOut   = fopen(FileOut,'w');
%
if TeXheader == 1
    fprintf(fidOut,'%s\n','\documentclass[a4paper,12pt]{report}');
    fprintf(fidOut,'%s\n','\usepackage{graphicx}');
    fprintf(fidOut,'%s\n','\usepackage{amsmath}');
    fprintf(fidOut,'%s\n','\usepackage{amsthm}');    
    fprintf(fidOut,'%s\n','\usepackage{amsfonts}');    
    fprintf(fidOut,'%s\n','\usepackage{longtable}');
    fprintf(fidOut,'%s\n','\begin{document}');
    str  = sprintf('%sitle{%s}','\t',Title);
    fprintf(fidOut,'%s\n',str);
    str  = sprintf('%suthor{%s}','\a',AuthorName);
    fprintf(fidOut,'%s\n',str);  
    fprintf(fidOut,'%s\n','\maketitle');
end
%
fprintf(fidOut,'%s\n','\begin{verbatim}');
%
tline    = fgetl(fidIn);
%
while ischar(tline)
    cont = cont + 1;
    if ~isempty(tline)
%         if tline(1) == '%'
%             fprintf(fidOut,'%s\n',tline(1:end));
%         end
        tmp1 = strfind(tline,'%c');
        if ~isempty(tmp1)
   %         disp(sprintf('Linha %d eh um comentario',cont))
            fprintf(fidOut,'%s\n',tline(tmp1+2:end));
        end
    end
      tline = fgetl(fidIn);
end
%
fprintf(fidOut,'%s\n','\end{verbatim}');
%
if TeXheader == 1
    fprintf(fidOut,'%s\n','\end{document}');
end

%
fclose(fidIn);
fclose(fidOut);
