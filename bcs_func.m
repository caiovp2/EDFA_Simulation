function res = bcs_func( ya,yb,Signal,Pump,ASE )

res (1) = ya(1)-Signal.Power; % Signal -> L(0) = Signal.Power
res (2) = ya(2)-Pump.Power;   % Pump   -> L(0) = Pump.Power
res (3) = ya(3)-ASE.Power;    % ASE+   -> L(0) = ASE.Power
res (4) = yb(4)-ASE.Power;    % ASE-   -> L(L) = ASE.Power

end

