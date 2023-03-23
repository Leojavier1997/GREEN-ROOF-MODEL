classdef Air
   properties (Access = public)
        
        Cp=1005.6;     % specific heat air
        R=286;      % gas constant of air [J/kg K]
        dens=40;  % molar density of air
        k=0.02554;    % aire [W/m K]
        gamma=0.06884;        
        Kv=0.4;
        Sch=0.63;      % turbulent Schmidt number
        Pr=0.71;      % turbulent Prandtl number

   end 
    
   methods
       function d = density(obj,Pa,Tair)
           d = Pa/(obj.R*Tair); 
       end
       
       function mu = kinematicViscosity(obj,T)
           % Function obtained from fitting data 
           % available in http://www.engineeringtoolbox.com/dry-air-properties-d_973.html           
           mu = 5.409413e-11*T*T+7.490993e-8*T-1.163350e-5;
       end
   end
end