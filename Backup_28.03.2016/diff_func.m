function dP = diff_func(z,P,options,Fiber,Signal,Pump,ASE,h,m,c)

A = [Signal.Absorption , Pump.Absorption, ASE.Absorption, ASE.Absorption];
G = [Signal.Gain , Pump.Gain, ASE.Gain, ASE.Gain];
v = c./[Signal.Wavelength , Pump.Wavelength, ASE.Wavelength, ASE.Wavelength];
u = [Signal.u , Pump.u, ASE.u, -ASE.u];
L = (log(10)/10).*[Signal.l , Pump.l, ASE.l, ASE.l];

s = length(Signal.Power);
p = length(Pump.Power);
a = 2*length(ASE.Power);
dP = zeros(s+p+a,1);
                         
mm = 0;
nn = 1;

for kk = 1 : s+p+a
 mm = P(kk,1)*A(1,kk)/(h*v(1,kk)*Fiber.zeta) + mm;
 nn = P(kk,1)*(A(1,kk)+G(1,kk))/(h*v(1,kk)*Fiber.zeta) + nn;
end

N = mm/nn;

%% TEST

% SIGNAL
 dP(1,1) = u(1,1)*(A(1,1)+G(1,1))*N*P(1,1)  - ...
           u(1,1)*(A(1,1)+L(1,1))*P(1,1)...
                  +u(1,1)*G(1,1)*N*m*h*v(1,1)*ASE.BW;
% PUMP
 dP(2,1) = u(1,2)*(A(1,2)+G(1,2))*N*P(2,1)  - ...
           u(1,2)*(A(1,2)+L(1,2))*P(2,1)...
                  +u(1,2)*G(1,2)*N*m*h*v(1,2)*ASE.BW;
% ASE+
 dP(3,1) = u(1,3)*(A(1,3)+G(1,3))*N*P(3,1)  + ...
           u(1,3)*G(1,3)*N*m*h*v(1,3)*ASE.BW - ...
           u(1,3)*(A(1,3)+L(1,3))*P(3,1);  
% ASE-
 dP(4,1) = u(1,4)*(A(1,4)+G(1,4))*N*P(4,1)  + ...
           u(1,4)*G(1,4)*N*m*h*v(1,4)*ASE.BW - ...
           u(1,4)*(A(1,4)+L(1,4))*P(4,1);  
 