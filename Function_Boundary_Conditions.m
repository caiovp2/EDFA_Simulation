function res = Function_Boundary_Conditions( ya,yb,Signal,Pump,ASE )

s = length(Signal.Wavelength);
p = length(Pump.Wavelength);
a = 2*length(ASE.Wavelength);
res = zeros(1,s+p+a);

res (1) = ya(1)-Signal.Power; % Signal -> L(0) = Signal.Power
res (2) = ya(2)-Pump.Power;   % Pump   -> L(0) = Pump.Power

for kk=s+p+1:s+p+a/2
res (kk) = ya(kk)-ASE.Power;    % ASE+   -> L(0) = ASE.Power
end

for kk=s+p+(a/2)+1:s+p+a
res (kk) = yb(kk)+ASE.Power;    % ASE-   -> L(L) = ASE.Power
end

end

