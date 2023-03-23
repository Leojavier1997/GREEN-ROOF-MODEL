classdef WeatherDataLine
    properties (Access = public)
        %% measured variables
        Tair % Air temperature
        RH % Relative Humidity
        U % Wind speed
        rainfall %        
        R_sh %Incomming solar Radiation
        %VWC % Volumetric water content in soil
        %Tin % Interior temperature (below the ceiling)
        Pa % Atmospheric pressure
        irrigation
        IR % Infrared radiation
        %% derived variables
        Tsky % Sky temperature
        
        
    end
    methods
        function obj = WeatherDataLine(data)
            
            %c = Constants;
            
            %% Measured
            obj.rainfall = data(1); % Rainfall             
            obj.R_sh = max(data(2),0); % Solar irradiance
            obj.Tair = data(3) + 273.15; % Exterior Air Temperature
            obj.RH = data(4)/100; % Exterior air RH 
            obj.U = data(5); % Wind speed 
            obj.Pa = 101300; %Pa... assumed constant                        
            obj.irrigation = data(6); %data(8)/10;
                       
            
            %% Derived
            
            
            consts = Constants();
            
            % Assume clear sky
            b=18.678;
            c = 275.14;
            d = 234.5;
            g = log(obj.RH*exp((b-(obj.Tair-273.15)/d)*((obj.Tair-273.15)/(c+(obj.Tair-273.15)))));            
            dewPointT = c*g/(b-g);
            skyEmissivity = 0.787+0.764*log((dewPointT+273.15)/274);
            obj.IR = skyEmissivity*consts.SB*obj.Tair^4;
            obj.Tsky = (obj.IR/consts.SB)^0.25;            
        end
    end
    
end