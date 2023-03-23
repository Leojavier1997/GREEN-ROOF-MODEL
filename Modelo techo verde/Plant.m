classdef Plant
   properties (Access = public)
        rho=.06;     %sw reflectivity plant
        em=0.98;    %emissivity planta
        ks=0.14; % Extintion coefficient   
        ks_ir = 0.14; % IR extintion coefficient
        LAI=4;        %leaf area index A FEB Tabares  
        
        k=0.5;      % thermal conductivity planta
        rsmin=500;      %[s/m] resistencia estomatica 
        height=0.55;   %  height of plant [m]
        z=2;       % [m] height of wind measurements
        
        Zog = 0.001; % VARIA SEGUN SMOOOTH
        
        %% Sailor requirements
        %fc=0.01;    % fractional vegetation coverage %
                
   end
   methods
        function obj = Plant()
            
            
        end       
               
       
        function s = fc(obj)
            s = 0.9 - 0.7*exp(-0.75*obj.LAI);
        end
        
        function r = tau_fsol(obj)
           r = exp(-obj.ks*obj.LAI);
        end
       
        function r = tau_fir(obj)
           r = exp(-obj.ks_ir*obj.LAI);
        end
        
        function r = foliage_rho(obj)
           r = (1-obj.tau_fsol)*obj.rho; 
        end
   end
   
end