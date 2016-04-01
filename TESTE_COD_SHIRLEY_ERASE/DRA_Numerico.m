%% Modelo numérico CW, prepara variáveis para rodar função bvp4c.
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
% Saída:
% resultados_num: Valores obtidos com modelo numérico
%                P_evol: [2*p+4*n x k] %evolução da potência (bombeio,bombeio esp, sinal, sinal esp, ASE, ASE esp)
%               lambda: [1xm] % lambada dos bombeios
%                   P : [1xm] % Potência dos bombeios em z = L 
%   Ganho_on_off_medio: double %Sem atenuação
%          Ganho_Medio: double %Com atenuação
%         Ganho_on_off: [1xn] %Sem atenuação
%               ripple: double %ripple com atenuação
%
% resultados_an: Valores obtidos com modelo analítico
%               lambda: [1xm] % lambada dos bombeios
%                   P : [1xm] % Potência dos bombeios em z = L 
%               P_evol: [1xm] % potência dos bombeios em z = 0
%   Ganho_on_off_medio: double %Sem atenuação
%          Ganho_Medio: double %Com atenuação
%         Ganho_on_off: [1xn] %Sem atenuação
%               ripple: double %ripple com atenuação
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [resultados_num,resultados_an] = DRA_Numerico(bombeio, sinal, fibra)
format long

%% Roda modelo analítico para determinar chute da Potência dos bombeios em z = 0
[ripple_an,Ganho_Medio_an,Ganho_on_off_medio_an,GA_sinaldB_an, Ppumpch] = DRA_Analitico_num(bombeio,sinal,fibra);

resultados_an.lambda = bombeio.lambda;
resultados_an.P = bombeio.P;
resultados_an.P_evol = Ppumpch;
resultados_an.ripple = ripple_an;
resultados_an.Ganho_Medio = Ganho_Medio_an;
resultados_an.Ganho_on_off_medio = Ganho_on_off_medio_an;
resultados_an.Ganho_on_off = GA_sinaldB_an;

if length(sinal.P) == 1
    sinal.P = sinal.P*ones(1,sinal.N_Sinais);
end

%% Preparação de variáveis para o método numérico
%% CHUTES
Ppumpbch=1e-4*ones(1,bombeio.N_Bombeios);   % Chutes dos pumps espalhados
Posch=sinal.P;       % Chutes dos sinais 
Pobch=1e-6*ones(1,sinal.N_Sinais);          % Chutes dos sinais espalhados
Pasech=zeros(1, sinal.N_Sinais);
Pasebch=zeros(1, sinal.N_Sinais);

%% Dados invariáveis
x = linspace(0,fibra.L,50);         % Vetor posição
Pos = sinal.P;

format long
%% Determina parâmetros do amplificador Raman
Param = Parametros(bombeio,sinal,fibra);

%% BVP - Inicialização dos chutes
solinit=bvpinit(x,[Ppumpch Ppumpbch Posch Pobch Pasech Pasebch]);  %Ppo Ppespl Pos Pbl Pase Paseb
options = bvpset('RelTol',1e-9,'AbsTol',1e-9,'Nmax',1e3);

%% Solução da Equação Diferencial Ordinária

num_ODE = @(x,y) numerico_ODE(x,y,bombeio,sinal,fibra,Param);
num_BC = @(ya,yb) numerico_BC(ya,yb,bombeio,sinal);

sol = bvp4c(num_ODE,num_BC,solinit,options);
%sol = bvp5c(num_ODE,num_BC,solinit,options); %Erro(matriz esparsa)

P = deval(sol,x);
sol.y;

%% Resultados

% Evolução da potência(bombeios, bombeios espalhado, sinal, sinal
% espalhado, ASE, ASE espalhada)
resultados_num.P_evol = P;

lx=length(x);

%% Evolução da potência dos sinais
for i=1:sinal.N_Sinais
    Psinais(i,:) = P(2*bombeio.N_Bombeios+i,:);                  
    Psinais_ASE(i,:) = P(2*bombeio.N_Bombeios+i,:)+P(2*(bombeio.N_Bombeios+sinal.N_Sinais)+i,:);
    Psinais_esp(i,:) = P(2*bombeio.N_Bombeios+sinal.N_Sinais+i,:);    
     
    % Net Gain com Sinal/Ruído
    NGsr(i,:)=(P(2*bombeio.N_Bombeios+i,:)+P(2*(bombeio.N_Bombeios+sinal.N_Sinais)+i,:))./Pos(i);  %evoluçao do NG com o ruido somado ao sinal
    NGfinalsr(i)=NGsr(i,lx);
    NGdBsr(i,:)=10*log10((P(2*bombeio.N_Bombeios+i,:)+P(2*(bombeio.N_Bombeios+sinal.N_Sinais)+i,:))./Pos(i));   %Net Gain dB considerando sinal mais ruido
    NGfinaldBsr(i)=NGdBsr(i,lx);

    % Relação Sinal/Ruído Óptica
    OSNRfinalsr(i)=10*log10(((P(2*bombeio.N_Bombeios+i,lx))+P(2*(bombeio.N_Bombeios+sinal.N_Sinais)+i,lx))/(P(2*(bombeio.N_Bombeios+sinal.N_Sinais)+i,lx))); %OSNR considerando o sinal na saida mais o ruido
    
    % Potência ASE
    PasedB(i)=10*log10(P(2*(bombeio.N_Bombeios+sinal.N_Sinais)+i,lx)/1e-3);
end

NGfimsr=NGfinaldBsr;
ripplesr=max(NGfimsr)-min(NGfimsr);
NGmediosr=mean(NGfimsr);
Ganho_Medio_OnOff = mean(NGfimsr + fibra.alfasdBkm*fibra.L/1e3);

resultados_num.lambda = bombeio.lambda;
resultados_num.P = bombeio.P;
resultados_num.Ganho_on_off_medio = Ganho_Medio_OnOff;
resultados_num.Ganho_Medio=NGmediosr;
resultados_num.Ganho_on_off = NGfimsr + fibra.alfasdBkm*fibra.L/1e3;
resultados_num.ripple=ripplesr;
