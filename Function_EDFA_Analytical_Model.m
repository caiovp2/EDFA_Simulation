function  [Ps,Pp]= Function_EDFA_Analytical_Model(Fiber,Signal,Pump,h,c)
 
 %Qin = P/hv = P*lambda/hc
 Qin = [Signal.Power*Signal.Wavelength Pump.Power*Pump.Wavelength]./(h*c);
 A = [Signal.Absorption Pump.Absorption];
 G = [Signal.Gain Pump.Gain];
 
 implicit = @(x,y) Function_EDFA_Analytical_Implicit_Equation(x,Fiber,Qin,A,G);
 xo = [0 sum(Qin)];
 Qout = fzero(implicit,xo);

 Ps = Qin(1,1)*exp( ( A(1,1)+G(1,1) )*( Qin(1,1)+Qin(1,2)-Qout )/Fiber.zeta-...
                               A(1,1)*Fiber.Length )*h*c/Signal.Wavelength;
                           
 Pp = Qin(1,2)*exp( ( A(1,2)+G(1,2) )*( Qin(1,1)+Qin(1,2)-Qout )/Fiber.zeta-...
                                 A(1,2)*Fiber.Length )*h*c/Pump.Wavelength;

end

