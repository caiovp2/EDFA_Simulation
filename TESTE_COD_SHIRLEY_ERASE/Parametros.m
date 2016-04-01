%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Entradas:
% bombeio =                  sinal =                  fibra =
%    N_Bombeios: int(m)          lambda: [1xn]                  L: double 
%        lambda: [1xm]                P: [1xn]              alfas: [1xn] 
%             P: [1xm]         N_Sinais: int(n)             alfap: [1xm] 
%           Bwp: double             Bws: double         alfasdBkm: [1xn]   
%           FPL: double                                 alfapdBkm: [1xm]  
%                                                     Crpicosinal: [1xn]  
%                                                      Crpicopump: [1xm]
%                                                        Crnormal: [1xK]
%                                                         sepfreq: [1xK]
%                                                           Aeffs: [1xn]  
%                                                           Aeffp: [1xm]        
%                                                              KR: double 
%                                                              NA: double 
%                                                              no: double 
%
% Saídas:
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

function Param = Parametros(bombeio,sinal,fibra)

np = bombeio.N_Bombeios;
ns = sinal.N_Sinais;

fp = (2.99792458e8./bombeio.lambda)*10^9;  %Hz
fs = (2.99792458e8./sinal.lambda)*10^9;  %Hz
wp=2*pi*fp;
ws=2*pi*fs;

%Parametros p/ calculo da dependencia da temperatura
hP= 6.626068e-34;     % constante de Planck em (Js)
Kb= 1.3806503e-23;    % constante de Boltzmann (J/K)
Tf= 300;              % temperatura da fibra (K)
c = 2.99792458e8;

%% Banda de Ruído de Referência Be (Eq. 2.27 e 2.28 - Tese da Shirley)
% Cálculo do ruído Rayleigh e ASE nas freq. dos SINAIS %
Bes = (c./sinal.lambda.^2).*sinal.Bws.*10^9; % largura de banda around fs (Hz)
Param.Nes= hP.*fs.*Bes ; % energia gerada pelo sinal devido à emissão espontânea   

% Cálculo do ruído Rayleigh e ASE nas freq. dos PUMPS 
% Só existe ASE de pump para multipumps%
Bep = (c./bombeio.lambda.^2).*bombeio.Bwp.*10^9; % largura de banda around fs (Hz)
Param.Nep= hP.*fp.*Bep;   % energia gerada pelo sinal devido à emissão espontânea   

%% Coeficiente de Rayleigh (Eq. 2.21 - Tese da Shirley)
%fibra.KR; %60;   % (dB/km) constante que depende do material dopante da fibra
%fibra.NA; %0.14; % abertura numerica da fibra
%fibra.no; %1.468;  % indice de refracao do nucleo da fibra
deltaindice = fibra.NA^2/(2*fibra.no);

% Perda de background do espalhamento Rayleigh
alfaRSdBkm = 0.63+fibra.KR*deltaindice;     %(dB/km)
alfaRS = 10^(alfaRSdBkm/10)*1e-3;     %(1/m)

Param.Rss = ((1000e-9)^4*alfaRS)./(4*pi*fibra.no^2.*(sinal.lambda*1e-9).^2.*fibra.Aeffs);  %1/m
Param.Rpp = ((1000e-9)^4*alfaRS)./(4*pi*fibra.no^2.*(bombeio.lambda.*1e-9).^2.*fibra.Aeffp);  %1/m

%% Calculo do ganho de Raman

pontos = length(fibra.Crnormal);

% Crnovo entre pumps e sinais
deltafps = zeros(np,ns);
Crnovops = zeros(np,ns);
Crps = zeros(np,ns);

for k=1:sinal.N_Sinais
    for m=1:bombeio.N_Bombeios
        deltafps(m,k)=abs(fp(m)-fs(k));
        for i=1:pontos-1
            if fibra.sepfreq(i) <= abs(deltafps(m,k)) && fibra.sepfreq(i+1)>=abs(deltafps(m,k))
                Crnovops(m,k)=fibra.Crnormal(i)+(((abs(deltafps(m,k)))-fibra.sepfreq(i))*...
                    (fibra.Crnormal(i+1)-fibra.Crnormal(i))/(fibra.sepfreq(i+1)-fibra.sepfreq(i)));
                break;
            end
        end
        Crps(m,k)=fibra.Crpicopump(m)*Crnovops(m,k); %(m/W)
    end
end

deltavps = zeros(np,ns);
Param.gsps = zeros(np,ns);
Param.gpps = zeros(np,ns);
Param.Taseps= zeros(np,ns);

for i=1:np
    for j=1:ns
        deltavps(i,j)=abs(fp(i)-fs(j));
        Param.gsps(i,j)=Crps(i,j)/bombeio.FPL;
        Param.gpps(i,j)=(wp(i)./ws(j)).*(Crps(i,j)/bombeio.FPL);
        %Distribuição de Bose-Einstein (Eq. 2.26)
        Param.Taseps(i,j)=1+(1/((exp(hP*deltavps(i,j)/(Kb*Tf)))-1)); 
    end
end

%Crnovo entre pumps
deltafpp = zeros(np);
Crnovopp = zeros(np);
Crpp = zeros(np);
for k=1:np
    for m=1:np
        deltafpp(m,k)=abs(fp(m)-fp(k)); 
         for i=1:pontos-1
             if fibra.sepfreq(i) <= abs(deltafpp(m,k)) && fibra.sepfreq(i+1) >= abs(deltafpp(m,k))
                Crnovopp(m,k)=fibra.Crnormal(i)+(((abs(deltafpp(m,k)))-fibra.sepfreq(i))*...
                    (fibra.Crnormal(i+1)-fibra.Crnormal(i))/(fibra.sepfreq(i+1)-fibra.sepfreq(i)));
                break;
            end
        end
       Crpp(m,k)=fibra.Crpicopump(m).*Crnovopp(m,k);
   end
end

Param.gspp = zeros(np);
Param.gppp = zeros(np);
Param.Tasepp = zeros(np);
for i=1:np
    for k=1:np
        if i == k
            deltafpp(i,k)=0;
            Param.gspp(i,k)=0;
            Param.gppp(i,k)=0;
            Param.Tasepp(i,k)=0;
        elseif i < k
            Param.gppp(i,k)=(wp(i)./wp(k)).*(Crpp(i,k)/bombeio.FPL);   
            Param.gspp(i,k)=0;
            Param.Tasepp(i,k)=1+(1/((exp(hP*deltafpp(i,k)/(Kb*Tf)))-1));
        else %i > k
            Param.gspp(i,k)=Crpp(i,k)/bombeio.FPL;        
            Param.gppp(i,k)=0;
            Param.Tasepp(i,k)=1+(1/((exp(hP*deltafpp(i,k)/(Kb*Tf)))-1));              
        end
    end
end
% Crnovo entre sinais
deltafss = zeros(ns);
Crnovoss = zeros(ns);
Crss = zeros(ns);
for k=1:ns
    for m=1:ns
        deltafss(m,k)=abs(fs(m)-fs(k));
         for i=1:pontos-1
             if fibra.sepfreq(i) <= abs(deltafss(m,k)) && fibra.sepfreq(i+1) >= abs(deltafss(m,k))
                Crnovoss(m,k)= fibra.Crnormal(i)+(((abs(deltafss(m,k))) - fibra.sepfreq(i))*...
                (fibra.Crnormal(i+1)-fibra.Crnormal(i))/(fibra.sepfreq(i+1)-fibra.sepfreq(i)));
                break;
            end
        end
       Crss(m,k)=fibra.Crpicosinal(m)*Crnovoss(m,k);
   end
end

Param.gsss = zeros(ns);
Param.gpss = zeros(ns);
Param.Tasess = zeros(ns);
for i=1:ns
    for k=1:ns
        if i == k 
            deltafss(i,k)=0;
            Param.gsss(i,k)=0;
            Param.gpss(i,k)=0;
            Param.Tasess(i,k)=0;
        elseif i < k
            Param.gpss(i,k)=(ws(i)./ws(k)).*(Crss(i,k)/bombeio.FPL); 
            Param.gsss(i,k)=0;
            Param.Tasess(i,k)=1+(1/((exp(hP*deltafss(i,k)/(Kb*Tf)))-1));       
        else %i > k
            Param.gsss(i,k)=Crss(i,k)/bombeio.FPL;           
            Param.gpss(i,k)=0;
            Param.Tasess(i,k)=1+(1/((exp(hP*deltafss(i,k)/(Kb*Tf)))-1));
        end
    end
end