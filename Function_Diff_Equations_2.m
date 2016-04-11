function dP = Function_Diff_Equations_2(z,P,Fiber,Signal,Pump,ASE,h,m,c)

minargs  = 9; 
maxargs  = 9;
narginchk(minargs, maxargs) % Checks if the number of inputs is between
                            % minargs and maxargs
                            

A = [Signal.Absorption  , Pump.Absorption , ASE.Absorption , ASE.Absorption ];
G = [Signal.Gain  , Pump.Gain , ASE.Gain , ASE.Gain ];
v = c./[Signal.Wavelength  , Pump.Wavelength , ASE.Wavelength , ASE.Wavelength ];
u = [Signal.u , Pump.u, ASE.u, -ASE.u];
L = (log(10)/10).*[Signal.l , Pump.l, ASE.l, ASE.l];

s = length(Signal.Wavelength);
p = length(Pump.Wavelength);
a = 2*length(ASE.Wavelength);
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
           u(1,1)*(A(1,1)+L(1,1))*P(1,1);
% PUMP
 dP(2,1) = u(1,2)*(A(1,2)+G(1,2))*N*P(2,1)  - ...
           u(1,2)*(A(1,2)+L(1,2))*P(2,1);
       
for kk=s+p+1:s+p+a/2       
% ASE+
 dP(kk,1) = u(1,kk)*(A(1,kk)+G(1,kk))*N*P(kk,1)  + ...
           u(1,kk)*G(1,kk)*N*m*h*v(1,kk)*ASE.BW - ...
           u(1,kk)*(A(1,kk)+L(1,kk))*P(kk,1);  
end

for kk=s+p+(a/2)+1:s+p+a
% ASE-
 dP(kk,1) = u(1,kk)*(A(1,kk)+G(1,kk))*N*P(kk,1)  + ...
           u(1,kk)*G(1,kk)*N*m*h*v(1,kk)*ASE.BW - ...
           u(1,kk)*(A(1,kk)+L(1,kk))*P(kk,1);
end
