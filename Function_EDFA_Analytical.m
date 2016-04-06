function y = Function_EDFA_Analytical(Qout,Fiber,Signal,Pump,h,c)

% Qin = P/hv = P*lambda/hc
Qin = [Signal.Power*Signal.Wavelength Pump.Power*Pump.Wavelength]./(h*c);
A = [Signal.Absorption Pump.Absorption];
G = [Signal.Gain Pump.Gain];
L = (log(10)/10).*[Signal.l Pump.l];

y = Qin(1,1)*exp( ( A(1,1)+G(1,1) )*( Qin(1,1)+Qin(1,2)-Qout )/Fiber.zeta -A(1,1)*Fiber.Length ) +...
    Qin(1,2)*exp( ( A(1,2)+G(1,2) )*( Qin(1,1)+Qin(1,2)-Qout )/Fiber.zeta -A(1,2)*Fiber.Length ) -...
    Qout;

end