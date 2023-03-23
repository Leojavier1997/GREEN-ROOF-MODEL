classdef Substrate    
    properties (Access = public)        
        rho=0.1853;     %sw reflectivity substrate
        depth=0.10;     % depth substrate  (ex t)
        Lsubm=0.05;    % depth substrate middepth sensor
        em = 0.96;
        dens =1100; %kg/m^3
        Cp=1000; % J/kgK;
        k=0.13; %W/mK;
        VWCsat=0.8;     % VWC saturated
        VWCfc = 0.35;
        VWCwilting=0.02;    % VWC wilting point
        VWCresidual = 0.02;
        phi=0.85;           
    end
    
    methods
       function k = thermalConductivity(obj,VWC)
           k = (obj.k + VWC*0.591);           
       end
       
       function rho = density(obj,VWC)
           rho = obj.dens + VWC*1000;           
       end
       
       function rcp = rhoCp(obj,VWC)
           %r = obj.density(VWC);     
           % Alexantri, Jones, 2006
           rcp = (1-obj.VWCsat)*obj.dens*obj.Cp + VWC * 1000*4182;           
       end
       
    end
end