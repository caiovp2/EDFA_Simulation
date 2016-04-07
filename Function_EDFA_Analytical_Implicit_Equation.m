function y = Function_EDFA_Analytical_Implicit_Equation(Qout,Fiber,Qin,A,G)

y = Qin(1,1)*exp( ( A(1,1)+G(1,1) )*( Qin(1,1)+Qin(1,2)-Qout )/Fiber.zeta -A(1,1)*Fiber.Length ) +...
    Qin(1,2)*exp( ( A(1,2)+G(1,2) )*( Qin(1,1)+Qin(1,2)-Qout )/Fiber.zeta -A(1,2)*Fiber.Length ) -...
    Qout;

end