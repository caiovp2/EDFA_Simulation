% Arquivo .m idêntico a versão compilada de numerico_BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Entradas:
% bombeio =                         sinal =                 
%    N_Bombeios: int(m)                 lambda: [1xn] (NU)  
%        lambda: [1xm] (NU)                  P: [1xn] (NU)  
%             P: [1xm] (NU)           N_Sinais: int(n)      
%           Bwp: double (NU)               Bws: double (NU) 
%           FPL: double (NU)                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = numerico_BC_m(ya,yb,bombeio,sinal)

Pp0 = 0;                % (W) Potencia dos pumps espalhados em z=0 (SEMPRE ZERO)
Pbl = 0;                % Potencia do sinal espalhado em L
Pasef0 = 0;             %potencia de ASE em z=0 (W)
PasebL = 0;             %potencia de ASE espalhada em z=L (W)


np = bombeio.N_Bombeios;
ns = sinal.N_Sinais; 
index_sinal_f           = 2*np;
index_sinal_b           = 2*np + ns;
index_ase_f             = 2*np + 2*ns;
index_ase_b             = 2*np + 3*ns;

Pump = yb(1:np) - bombeio.P';
Pumpespal = ya(np+1:index_sinal_f)-Pp0;
Sinal = ya(index_sinal_f+1:index_sinal_b) - sinal.P';
Sinalmenos = yb(index_sinal_b+1:index_ase_f) - Pbl; %corrigido, era 'ya'
ASE = ya(index_ase_f+1:index_ase_b) - Pasef0;
ASEespal = yb(index_ase_b+1:index_ase_b + ns) - PasebL;

res= [Pump; Pumpespal; Sinal; Sinalmenos; ASE; ASEespal];