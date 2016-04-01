% Arquivo .m idêntico a versão compilada de numerico_ODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Entradas:
% bombeio =                  sinal =                  fibra =
%    N_Bombeios: int(m)          lambda: [1xn] (NU)             L: double (NU)
%        lambda: [1xm] (NU)           P: [1xn] (NU)         alfas: [1xn] 
%             P: [1xm] (NU)    N_Sinais: int(n)             alfap: [1xm] 
%           Bwp: double (NU)        Bws: double (NU)    alfasdBkm: [1xn]  (NU) 
%           FPL: double (NU)                            alfapdBkm: [1xm]  (NU)
%                                                     Crpicosinal: [1xn]  (NU)
%                                                      Crpicopump: [1xm]  (NU)
%                                                        Crnormal: [1xK]  (NU)
%                                                         sepfreq: [1xK]  (NU)
%                                                           Aeffs: [1xn]  (NU)
%                                                           Aeffp: [1xm]  (NU)      
%                                                              KR: double (NU)
%                                                              NA: double (NU)
%                                                              no: double (NU)
%     Param = 
%         Nes: [1xn]
%         Nep: [1xm]
%         Rss: [1xn]
%         Rpp: [1xm]
%        gsps: [mxn]
%        gpps: [mxn]
%      Taseps: [mxn]
%        gspp: [mxm]
%        gppp: [mxm]
%      Tasepp: [mxm]
%        gsss: [nxn]
%        gpss: [nxn]
%      Tasess: [nxn]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydx = numerico_ODE_m(~,y,bombeio,sinal,fibra,Param)
np = bombeio.N_Bombeios;
ns = sinal.N_Sinais;

index_sinal_f           = 2*np;
index_sinal_b           = 2*np + ns;
index_ase_f             = 2*np + 2*ns;
index_ase_b             = 2*np + 3*ns;


%% Monta equações diferenciais
 alfas=fibra.alfas;
 alfap=fibra.alfap;
 Rs = Param.Rss;
 Rp = Param.Rpp;
 
 Pasef = zeros(1,ns);
 Paseb = zeros(1,ns);
 
 Psf = zeros(1,ns);
 Psb = zeros(1,ns);
 
for j=1:ns
    Nf1=0; Nf2=0; Nf3=0;
    Nase1=0;Nase2=0;
    Case1=0; Case2=0; 

    Nb1=0; Nb2=0; Nb3=0;
    Naseb1=0; Naseb2=0;
    Caseb1=0; Caseb2=0; 

    for k=1:np        
        %Ganho do sinal devido aos pumps +z e -z
        Nf1=Nf1+Param.gsps(k,j)*(y(k)+y(k+np))*y(index_sinal_f + j); 
        Nb1=Nb1+Param.gsps(k,j)*(y(k)+y(k+np))*y(index_sinal_b + j);
        %Ganho da ASE devido aos bombeios +z e -z   
        Nase1=Nase1+Param.gsps(k,j)*(y(k)+y(k+np))*y(index_ase_f + j); 
        Naseb1=Naseb1+Param.gsps(k,j)*(y(k)+y(k+np))*y(index_ase_b + j);
        %Geração da ASE em torno dos bombeios +z e -z
        Case1=Case1+Param.gsps(k,j)*(y(k)+y(k+np))*2*Param.Nes(j)*Param.Taseps(k,j); 
        Caseb1=Caseb1+Param.gsps(k,j)*(y(k)+y(k+np))*2*Param.Nes(j)*Param.Taseps(k,j); 
    end

    for m=1:ns  
        % Ganho do sinal devido aos outros sinais +z e -z
        Nf2=Nf2+Param.gsss(j,m)*(y(m+index_sinal_f)+y(m+index_sinal_b))*y(index_sinal_f + j); 
        Nb2=Nb2+Param.gsss(j,m)*(y(m+index_sinal_f)+y(m+index_sinal_b))*y(index_sinal_b + j);
        % Deplecao do sinal devido aos outros sinais
        Nf3=Nf3+Param.gpss(j,m)*((y(m+index_sinal_f)+y(m+index_sinal_b) + ...
            4*Param.Nes(m)*Param.Tasess(j,m))*y(index_sinal_f + j)); 
        Nb3=Nb3+Param.gpss(j,m)*((y(m+index_sinal_f)+y(m+index_sinal_b) + ...
            4*Param.Nes(m)*Param.Tasess(j,m))*y(index_sinal_b + j));
        %Ganho da ASE devido aos bombeios +z e -z   
        Nase2=Nase2+Param.gsss(j,m)*(y(m+index_sinal_f)+y(m+index_sinal_b))*y(index_ase_f + j);    
        Naseb2=Naseb2+Param.gsss(j,m)*(y(m+index_sinal_f)+y(m+index_sinal_b))*y(index_ase_b + j);
        %Geração da ASE em torno dos bombeios +z e -z
        Case2=Case2+Param.gsss(j,m)*(y(m+index_sinal_f)+y(m+index_sinal_b))*...
            2*Param.Nes(j)*Param.Tasess(j,m);
        Caseb2=Caseb2+Param.gsss(j,m)*(y(m+index_sinal_f)+y(m+index_sinal_b))*...
            2*Param.Nes(j)*Param.Tasess(j,m);
    end

    %% Eq. 2.31 (Tese da Shirley) - Para sinais
    % -Atenuação na fibra + Espalhamento duplo de Rayleigh + Ganho devido a 
    % bombeios + Ganho devido a sinais - Depleção devido aos sinais
    Psf(j) = -alfas(j)*y(index_sinal_f + j) + Rs(j)*y(index_sinal_b + j) + Nf1 + Nf2 - Nf3; % +z
    % Atenuação na fibra - Espalhamento duplo de Rayleigh - Ganho devido a 
    % bombeios - Ganho devido a sinais + Depleção devido aos sinais 
    Psb(j) = alfas(j)*y(index_sinal_b + j)- Rs(j)*y(index_sinal_f + j) - Nb1 - Nb2 + Nb3; % -z

    %% Eq. 2.31 (Tese da Shirley) - Para ASE 
    % -Atenuaçao da fibra + Espalhamento duplo de Rayleigh + Ganho da ASE 
    % devido aos bombeios + Ganho da ASE devido aos sinais + Amplificação
    % da ASE em torno dos bombeios + Amplificação da ASE em torno dos
    % Sinais
    Pasef(j)=-alfas(j)*y(index_ase_f + j) + Rs(j)*y(index_ase_b + j) + Nase1 + Nase2 + Case1 + Case2; %+z
    % +Atenuaçao da fibra - Espalhamento duplo de Rayleigh - Ganho da ASE 
    % devido aos bombeios - Ganho da ASE devido aos sinais - Amplificação
    % da ASE em torno dos bombeios - Amplificação da ASE em torno dos
    % Sinais
    Paseb(j)= +alfas(j)*y(index_ase_b + j) - Rs(j)*y(index_ase_f + j) - Naseb1 - Naseb2 - Caseb1 - Caseb2; %-z
end
 
Ppf = zeros(1,np);
Ppb = zeros(1,np);

for k=1:np
    Np1f = 0; Np2f = 0; Np3f = 0;
    Np1b = 0; Np2b = 0; Np3b = 0;

    for j=1:ns
        %Depleção dos bombeios devido aos sinais(+z e -z), ASE(+z e -z) e
        %Depleção devido a amplificação da ASE nas menores frequências(+z e -z)
        Np1f = Np1f + Param.gpps(k,j)*(y(index_sinal_f + j)+y(index_sinal_b + j)+y(index_ase_f + j)+y(index_ase_b + j)+...
            4*Param.Nes(j)*Param.Taseps(k,j))*y(k); 
        Np1b = Np1b + Param.gpps(k,j)*(y(index_sinal_f + j)+y(index_sinal_b + j)+y(index_ase_f + j)+y(index_ase_b + j)+...
            4*Param.Nes(j)*Param.Taseps(k,j))*y(k+np);
    end

    for n=1:np
        %Deplecao dos bombeios devido aos outro bombeios(+z e -z) e 
        %Depleção devido a amplificação da ASE nas menores frequências(+z e -z)
        Np2f = Np2f + Param.gppp(k,n)*(y(n)+y(n+np)+ 4*Param.Nep(n)*Param.Tasepp(k,n))*y(k); 
        Np2b = Np2b + Param.gppp(k,n)*(y(n)+y(n+np)+ 4*Param.Nep(n)*Param.Tasepp(k,n))*y(k+np);
        
        %Ganho dos bombeios devidos aos outros bombeios (+z  e -z)
        Np3f = Np3f + Param.gspp(k,n)*(y(n)+y(n+np))*y(k); 
        Np3b = Np3b + Param.gspp(k,n)*(y(n)+y(n+np))*y(k+np); 
    end
    %% Eq. 2.31 (Tese da Shirley) - Para Bombeios 
    % Os bombeios são contra-propagantes
    % Atenuaçao da fibra - Espalhamento duplo de Rayleigh + Depleção dos
    % bombeios devidos aos outro sinais + Depleção devido aos outros
    % bombeios - Ganho devido aos outros bombeios
    Ppf(k) = alfap(k)*y(k) - Rp(k)*y(k+np) + Np1f + Np2f - Np3f; %-z
    % -Atenuaçao da fibra + Espalhamento duplo de Rayleigh - Depleção dos
    % bombeios devidos aos outro sinais - Depleção devido aos outros
    % bombeios + Ganho devido aos outros bombeios
    Ppb(k) = -alfap(k)*y(k+np) + Rp(k)*y(k) - Np1b - Np2b + Np3b; %+z
end

dydx = [Ppf'; Ppb'; Psf'; Psb'; Pasef'; Paseb'];