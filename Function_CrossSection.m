function [sigmaA , sigmaE] = Function_CrossSection(CrossSection,Laser,kk)

minargs  = 2; 
maxargs  = 3;
narginchk(minargs, maxargs) % Checks if the number of inputs is between
                            % minargs and maxargs
if(nargin<maxargs)
    kk=1;
end

sigmaA   = interp1(CrossSection.Wavelength,...
                   CrossSection.Absorption,...
                   Laser.Wavelength(1,kk),'nearest');
               
sigmaE   = interp1(CrossSection.Wavelength,...
                   CrossSection.Emission,...
                   Laser.Wavelength(1,kk),'nearest');               

end

