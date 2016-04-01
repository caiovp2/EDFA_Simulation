function dX = testando(t,X,options,Fiber,Signal,Pump,ASE,h,m)

dX = zeros(2,1);

dX(1,1) = 2*X(1,1)-2*X(2,1);
dX(2,1) = -3*X(1,1)+X(2,1);

% TESTANDO A FUNCAO ODE45
% EQUACAO DIFF
% X'=2X-2Y
% Y'=-3X+Y
%
%SOLUCAO
%X(T)=2EXP(-T)+3EXP(4T)
%Y(T)=3EXP(-T)-3EXP(4T)