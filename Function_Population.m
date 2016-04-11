function N2 = Function_Population(Fiber,Signal,Pump,ASE,h,c,vector)

% N2 =((h*Fiber.zeta)^-1)*(vector(1,:)*Signal.Absorption/(c/Signal.Wavelength)+...
%                          vector(2,:)*Pump.Absorption/(c/Pump.Wavelength)+...
%                          vector(3,:)*ASE.Absorption/(c/ASE.Wavelength)+...
%                          vector(4,:)*ASE.Absorption/(c/ASE.Wavelength))./...
%       (1 +((h*Fiber.zeta)^-1)*...
%      (vector(1,:)*(Signal.Absorption+Signal.Gain)/(c/Signal.Wavelength)+...
%       vector(2,:)*(Pump.Absorption+Pump.Gain)/(c/Pump.Wavelength)+...
%       vector(3,:)*(ASE.Absorption+ASE.Gain)/(c/ASE.Wavelength)+...
%       vector(4,:)*(ASE.Absorption+ASE.Gain)/(c/ASE.Wavelength)));
  
A = [Signal.Absorption  , Pump.Absorption , ASE.Absorption , ASE.Absorption ];
G = [Signal.Gain  , Pump.Gain , ASE.Gain , ASE.Gain ];
v = c./[Signal.Wavelength  , Pump.Wavelength , ASE.Wavelength , ASE.Wavelength ];

s = length(Signal.Wavelength);
p = length(Pump.Wavelength);
a = 2*length(ASE.Wavelength);

N2 = zeros(1,length(vector));

for aa=1:length(vector)
mm = 0;
nn = 1;

for kk = 1 : s+p+a
 mm = vector(kk,aa)*A(1,kk)/(h*v(1,kk)*Fiber.zeta) + mm;
 nn = vector(kk,aa)*(A(1,kk)+G(1,kk))/(h*v(1,kk)*Fiber.zeta) + nn;
end

N2(1,aa) = mm/nn;
end
end

