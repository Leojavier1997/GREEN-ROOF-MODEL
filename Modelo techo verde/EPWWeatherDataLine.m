classdef EPWWeatherDataLine
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
        IR     % Infrared radiation           
        Tsky % Sky temperature
        
        
    end
    methods
        function obj = EPWWeatherDataLine(data,Z,SiteWindBLHeight,SiteWindExp)
            
             WeatherFileWindModCoeff = 1.5863;
             
             % SKIPPED
             
             % --. N1, \field Year
             % --. N2, \field Month
             % --. N3, \field Day
             % --. N4, \field Hour
             % --. N5, \field Minute
             % --. A1, \field Data Source and Uncertainty Flags
             
            % 1. N6, \field Dry Bulb Temperature 
            obj.Tair = data(1)+273.15;
            % 2. N7, \field Dew Point Temperature                   
            % 3. N8, \field Relative Humidity  
            obj.RH = data(3)/100;
            
            obj.Pa = data(4);% 4. N9, \field Atmospheric Station Pressure
            % 5. N10, \field Extraterrestrial Horizontal Radiation
            % 6. N11, \field Extraterrestrial Direct Normal Radiation
            obj.IR = data(7);% 7. N12, \field Horizontal Infrared Radiation Intensity            
            obj.R_sh = data(8); % 8. N13, \field Global Horizontal Radiation            
            % 9. N14, \field Direct Normal Radiation
            % 10. N15, \field Diffuse Horizontal Radiation
            % 11. N16, \field Global Horizontal Illuminance
            % 12. N17, \field Direct Normal Illuminance
            % 13. N18, \field Diffuse Horizontal Illuminance                
            % 14. N19, \field Zenith Luminance
            % 15. N20, \field Wind Direction
            obj.U = data(16)*WeatherFileWindModCoeff*(Z/SiteWindBLHeight)^SiteWindExp;% 16. N21, \field Wind Speed
            % 17. N22, \field Total Sky Cover
            % 18. N23, \field Opaque Sky Cover (used if Horizontal IR Intensity missing)
            % 19. N24, \field Visibility
            % 20. N25, \field Ceiling Height
            % 21. N26, \field Present Weather Observation
            % 22. N27, \field Present Weather Codes
            % 23. N28, \field Precipitable Water
            % 24. N29, \field Aerosol Optical Depth
            % 25. N30, \field Snow Depth
            % 26. N31, \field Days Since Last Snowfall
            % 27. N32, \field Albedo
            obj.rainfall = data(28);% 28. N33, \field Liquid Precipitation Depth
            % 29. N34 \field Liquid Precipitation Quantity
            
            %% Not sure                      
            obj.irrigation = 0; 
            consts = Constants();
            obj.Tsky = (obj.IR/consts.SB)^0.25;
        end
    end
    
end