function dP = Function_Diff_Equations_1(z,P,Fiber,Signal,Pump,ASE,h,m,c)

A = [Signal.sigmaA  , Pump.sigmaA , ASE.sigmaA , ASE.sigmaA ];
E = [Signal.sigmaE  , Pump.sigmaE , ASE.sigmaE , ASE.sigmaE ];
O = [Signal.Overlap , Pump.Overlap , ASE.Overlap , ASE.Overlap ];
v = c./[Signal.Wavelength  , Pump.Wavelength , ASE.Wavelength , ASE.Wavelength ];
u = [Signal.u , Pump.u, ASE.u, -ASE.u];
L = (log(10)/10).*[Signal.l , Pump.l, ASE.l, ASE.l];

s = length(Signal.Power);
p = length(Pump.Power);
a = 2*length(ASE.Power);
dP = zeros(s+p+a,1);
                         
mm = 0;
nn = 1;

for kk = 1 : s+p+a
 mm = P(kk,1)*A(1,kk)*O(1,kk)*Fiber.tau/(h*v(1,kk)*Fiber.Aeff) + mm;
 nn = P(kk,1)*(A(1,kk)+E(1,kk))*O(1,kk)*Fiber.tau/(h*v(1,kk)*Fiber.Aeff) + nn;
end

N2 = (mm/nn)*Fiber.nt;
N1 = Fiber.nt - N2;

%% TEST

% SIGNAL
 dP(1,1) = u(1,1)*(N2*E(1,1)-N1*A(1,1))*O(1,1)*P(1,1)  - ...
           u(1,1)*L(1,1)*P(1,1);
% PUMP
 dP(2,1) = u(1,2)*(N2*E(1,2)-N1*A(1,2))*O(1,2)*P(2,1)  - ...
           u(1,2)*L(1,2)*P(2,1);
% ASE+
 dP(3,1) = u(1,3)*(N2*E(1,3)-N1*A(1,3))*O(1,3)*P(3,1)  + ...
           u(1,3)*E(1,3)*N2*O(1,3)*m*h*v(1,3)*ASE.BW - ...
           u(1,3)*L(1,3)*P(3,1);
% ASE-
 dP(4,1) = u(1,4)*(N2*E(1,4)-N1*A(1,4))*O(1,4)*P(4,1)  + ...
           u(1,4)*E(1,4)*N2*O(1,4)*m*h*v(1,4)*ASE.BW - ...
           u(1,4)*L(1,3)*P(4,1);
 
