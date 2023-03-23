% THIS IS BASED ON ENERGY PLUS CODE

classdef SailorSteady
    properties
        %% Current state
        T_plants
        T_substrate
        T_interior
        VWC % volumetric water content
        dt = 60; %seconds between measurements
        interface_heat_flux % Heat flux at the interface
        evaporation
        transpiration
        substrate_convection
        plant_convection
        
        et_mm_hour        
        
        %% Materials
        roof
        plant
        sub
        
        %% other                
        L1
    end
    
    methods
        function obj = SailorSteady(plant,sub,roof,n_layers,Area,dt)
            
            obj.roof = roof;
            obj.plant = plant;
            obj.sub = sub;
            
            %% Initialize
            obj.T_plants = 300; % K
            obj.T_substrate = 300; % K
            obj.T_interior = 300;
            
            %% Other                       
            obj.VWC = sub.VWC;
        end
        
        function obj = update(obj,data)
            
            % Load constants
            consts = Constants;
            air = Air;
            
            %% INPUTS
            R_sh=   data.R_sh;  % incoming SW radiation global horizontal
            Tair=   data.Tair;
            Tsky=   data.Tsky;
            RH=     data.RH;   % Relative humidity [%]
            %obj.T_interior =   data.T_interior;   % interior surface temperature
            Pa=     data.Pa;   % atmospheric pressure [Pa]
            U=      max(3,data.U);       % Wind speed [m/s]
            W=      obj.VWC;    %VWC            
            Latm = data.IR;   
            
            TT = [obj.T_plants, obj.T_substrate];
           
            Tf = TT(1);
            Tg = TT(2);
            
            % Basic definitions
            sigmaf = obj.plant.fc;
            e0 = 2; % This is constant on EPlus... not sure why
            rch = air.Sch;
            rche = air.Pr;
            
            e_1 = obj.plant.em + obj.sub.em - obj.sub.em*obj.plant.em;
            eair = RH*611.2*exp(17.67*(Tair-273.15)/(Tair-29.65));
            qa = 0.622*eair/(Pa-eair);
            %Rhoa = Pa/air.R*Tair;
            Rhoa = air.density(Pa,Tair);
            
            %% Energy budget on foliage layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Sensible heat in foliage layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eq. 4 - Sailor 2008
            Taf = (1-sigmaf)*Tair + sigmaf*(0.3*Tair+0.6*Tf+0.1*Tg);
            
            
            %Rhof = Pa/(air.R*Taf);
            Rhof = air.density(Pa,Tf);
            
            % Eq. 3 - Sailor 2008
            Rhoaf = (Rhoa+Rhof)/2;
            
            % Eq. 7 - Sailor 2008
            Zd = 0.701*obj.plant.height^0.979;
            
            % Eq. 8 - Sailor 2008... Limiting it to 0.02 was copied from
            % EPlus. It is not in the paper
            Zo = max(0.131*obj.plant.height^0.997,0.02);
            
            % Eq. 6 - Sailor 2008 (This term is called Cfhn in EPlus code)
            Chnf = (air.Kv/log((obj.plant.z-Zd)/Zo))^2;
            
            % Eq. 5 - Sailor 2008
            Waf = 0.83*sigmaf*U*sqrt(Chnf)+(1-sigmaf)*U;
            
            % Eq. 9 - Sailor 2008
            Cf = 0.01*(1+0.3/Waf);
            %Cf =   10*(1+0.3/Waf); %    bulk heat transfer coefficient

            % Eq. 2 - Sailor 2008... corrected by e0, as in EPlus code
            Hf = (e0+1.1*obj.plant.LAI*Rhoaf*air.Cp*Cf*Waf)*(Taf-Tf);
            
            obj.plant_convection = Hf;
            
            % Latent heat flux in the foliage layer           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Eq. 12 - Sailor 2008
            ra = 1/(Cf*Waf);
            
            % Eq. 11a - Sailor 2008
            f1 = 1/min(1,(0.004*R_sh+0.005)/(0.81*(0.004*R_sh+1.0)));
            
            % Unsure about these definitions... made to fit EPlus with
            MeanRootMoisture = W;
            MoistureResidual = obj.sub.VWCresidual;
            MoistureMax = obj.sub.VWCsat;
            
            % Eq. 11b - Sailor 2008... f2 does not become zero. This idea
            % comes from EnergyPlus code            
            if(MoistureMax == MoistureResidual)
                f2 = 1e-10;
            else
                f2 = 1/((MeanRootMoisture-MoistureResidual)/(MoistureMax-MoistureResidual));                
            end
            
            % Eq. 11c - Sailor 2008
            f3 = 1; %1/exp(-0*(esf-eair));
            
            % Eq. 10 - Sailor 2008
            r_s = obj.plant.rsmin*f1*f2*f3/obj.plant.LAI;
            
            % Eq. 13 - Sailor 2008
            rn = ra/(ra+r_s);
            
            % Described... not shown in equation in Sailor 2008                              
            % in EPlus code : Mg = Moisture/MoistureMax;
            Mg=     W/obj.sub.VWCsat;  
                       
            % Eq. 16 - Sailor 2008... the second part is from EPlus code
            Lef = 1.91846e6*(Tf/(Tf-33.91))^2;
            if(obj.T_plants < 273.15) % Less than 0 C
                Lef = 2.838e6;
            end
            
            % Saturation pressure... from EPlus.. aparently, common formula.
            esf = 611.2*exp(17.67*(Tf-273.15)/(Tf-29.65));            
            qsf = 0.622*esf/(Pa-esf);
            esg = 611.2*exp(17.67*((Tg-273.15)/(Tg-29.65))); %Pa saturation vapor pressure
			qsg = 0.622*esg/(Pa-esg); %Saturation mixing ratio at ground surface temperature.
            
            
            % Eq. 15 - Sailor 2008
            qaf =((1-sigmaf)*qa+sigmaf*(0.3*qa+0.6*qsf*rn+0.1*qsg*Mg))/(1-sigmaf*(0.6*(1.0-rn)+0.1*(1.0-Mg)));
			
            
            % Eq. 14 - Sailor 2008
            Lf = Lef*obj.plant.LAI*Rhoaf*Cf*Waf*rn*(qaf-qsf);
			
            %% Soil energy budget            
            %%%%%%%%%%%%%%%%%%%%%
             
                                    
            % Derivative of Saturation vapor pressure, which is used in the
            % calculation of derivative saturation specific humidity
            %Desf = 611.2*exp(17.67*((Tf-273.15)/(Tf-29.65)))*(17.67*(Tf-273.15)*(-1.0)*(Tf- 29.65)^(-2)+17.67/(Tf-29.65));
			%dqf = ( ( 0.622 * Pa ) / ( Pa - esf )^2 ) * Desf; %Derivative of saturation specific humidity
			
            %Latent heat vaporiobj.plant.ztion  at the ground temperature
			Leg = 1.91846e6*(Tg/(Tg-33.91))^2;
			%Check to see if ice is sublimating or frost is forming.
            if ( obj.T_substrate < 273.15 ) 
                Leg = 2.838e6; %per FASST documentation p.15 after eqn. 37.
            end
            
            %Desg = 611.2*exp(17.67*((Tg-273.15)/(Tg-29.65)))*(17.67*(Tg-273.15)*(-1.0)*(Tg- 29.65)^(-2) + 17.67/(Tg-29.65));
			%dqg = (0.622*Pa/(Pa-esg)^2)*Desg;
            
            %Final Ground Atmosphere Energy Balance
			%Density of air at the soil surface temperature
			%Rhog = Pa/(air.R*Tg);
            Rhog = air.density(Pa,Tg);
            
            % Eq. 19 - Sailor 2008
            Rhoag=(Rhoa+Rhog)/2;
            
            % Eq. 23 - Sailor 2008
			Rib=2*consts.g*obj.plant.z*(Taf-Tg)/((Taf+Tg)*Waf*Waf); %Richardson Number

			% Eq. 22 - Sailor 2008... as in the EPlus code (they are not
			% the same)
            if ( Rib < 0.0 ) 
				Gammah = ( 1.0 - 16.0 * Rib)^(-0.5);
            else 
                if ( Rib >= 0.19 ) 
					Rib = 0.19;
                end
				Gammah = ( 1.0 - 5.0 * Rib)^(-0.5 );
            end
            
            % Eq. 21 - Sailor 2008
            Chng = (air.Kv/log(obj.plant.z/obj.plant.Zog))^2 / rch; % bulk transfer coefficient near ground
			
            % Eq. 20 - Sailor 2008
            Chg = Gammah*((1-sigmaf)*Chng+sigmaf*Chnf);
			
            % Eq. 18 - Sailor 2008
            sheatg = e0 + Rhoag * air.Cp * Chg * Waf; % added the e0 windless correction
			Hg = sheatg * ( Taf - Tg ); % sensible flux TO soil (W/m^2) DJS Jan 2011 (eqn. 32 in Frankenstein 2004)
            
            % Eq. 20 - Sailor 2008
			Chne = ( air.Kv / log( obj.plant.z / obj.plant.Zog ) )^2 / rche;
			
            % Latent heat flux in soil layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Eq. 26 - Sailor 2008
            Ce = Gammah*((1-sigmaf)*Chne + sigmaf * Chnf); % this is in fact Ceg in eq (28)

			% Eq. 25 - Sailor 2008
            qg = Mg * qsg + ( 1.0 - Mg ) * qaf; 
            
            % Eq. 24 - Sailor 2008... with an extra Mg?
            Lg = Ce * Leg * Waf * Rhoag * ( qaf - qg ) * Mg; % In the FASST documentation there is NO Mg. However, in looking
			
            % Conduction
            r_roof = obj.roof.depth/obj.roof.k;
            r_substrate = obj.sub.depth/obj.sub.thermalConductivity(W);
            Qcond=(TT(2)-obj.T_interior)/(r_roof + r_substrate + consts.rsi_roof);
            
            
            
            %% UPDATE            
            obj.evaporation = Lg;
            obj.substrate_convection = Hg;
            obj.transpiration = Lf; % CORREGIR
            obj.plant_convection = Hf;
            
            obj.et_mm_hour = -(Lg/Leg + Lf/Lef)*3600;    
            
            %% Update interface_heat_flux
            obj.interface_heat_flux = Qcond;
            
            
        end
        
        function [Res] = ResFUN(obj,TT,data)
            
            % Load constants
            consts = Constants;
            air = Air;
            
            %% INPUTS
            R_sh=   data.R_sh;  % incoming SW radiation global horizontal
            Tair=   data.Tair;
            Tsky=   data.Tsky;
            RH=     data.RH;   % Relative humidity [%]
            %obj.T_interior =   data.T_interior;   % interior surface temperature
            Pa=     data.Pa;   % atmospheric pressure [Pa]
            U=      max(3,data.U);       % Wind speed [m/s]
            W=      obj.VWC;    %VWC            
            Latm = data.IR;   
            Tf = TT(1);
            Tg = TT(2);
            
            % Basic definitions
            sigmaf = obj.plant.fc;
            e0 = 2; % This is on EPlus... not sure why
            rch = air.Sch;
            rche = air.Pr;
            
            e_1 = obj.plant.em + obj.sub.em - obj.sub.em*obj.plant.em;
            eair = RH*611.2*exp(17.67*(Tair-273.15)/(Tair-29.65));
            qa = 0.622*eair/(Pa-eair);
            %Rhoa = Pa/air.R*Tair;
            Rhoa = air.density(Pa,Tair);
            
            %% Energy budget on foliage layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Sensible heat in foliage layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eq. 4 - Sailor 2008
            Taf = (1-sigmaf)*Tair + sigmaf*(0.3*Tair+0.6*Tf+0.1*Tg);
            
            %Rhof = Pa/(air.R*Taf);
            Rhof = air.density(Pa,Taf);
            
            % Eq. 3 - Sailor 2008
            Rhoaf = (Rhoa+Rhof)/2;
            
            % Eq. 7 - Sailor 2008
            Zd = 0.701*obj.plant.height^0.979;
            
            % Eq. 8 - Sailor 2008... Limiting it to 0.02 was copied from
            % EPlus. It is not in the paper
            Zo = max(0.131*obj.plant.height^0.997,0.02);
            
            % Eq. 6 - Sailor 2008 (This term is called Cfhn in EPlus code)
            Chnf = (air.Kv/log((obj.plant.z-Zd)/Zo))^2;
            
            % Eq. 5 - Sailor 2008
            Waf = 0.83*sigmaf*U*sqrt(Chnf)+(1-sigmaf)*U;
            
            % Eq. 9 - Sailor 2008
            Cf = 0.01*(1+0.3/Waf);            
            
            % Eq. 2 - Sailor 2008... corrected by e0, as in EPlus code
            Hf = (e0+1.1*obj.plant.LAI*Rhoaf*air.Cp*Cf*Waf)*(Taf-Tf);
            
            % Latent heat flux in the foliage layer           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Eq. 12 - Sailor 2008
            ra = 1/(Cf*Waf);
            
            % Eq. 11a - Sailor 2008
            f1 = 1/min(1,(0.004*R_sh+0.005)/(0.81*(0.004*R_sh+1.0)));
            
            % Unsure about these definitions... made to fit EPlus with
            MeanRootMoisture = W;
            MoistureResidual = obj.sub.VWCresidual;
            MoistureMax = obj.sub.VWCsat;
            
            % Eq. 11b - Sailor 2008... f2 does not become zero. This idea
            % comes from EnergyPlus code            
            if(MoistureMax <= MoistureResidual)
                f2 = 1e-10; % COMPARAR CON MAIN REPORT
            else
                f2 = 1/((MeanRootMoisture-MoistureResidual)/(MoistureMax-MoistureResidual));                
            end
            
            % Eq. 11c - Sailor 2008
            f3 = 1; %1/exp(0*(esf-air));
            
            % Eq. 10 - Sailor 2008
            r_s = obj.plant.rsmin*f1*f2*f3/obj.plant.LAI;
            
            % Eq. 13 - Sailor 2008
            rn = ra/(ra+r_s);
            
            % Described... not shown in equation in Sailor 2008                              
            % in EPlus code : Mg = Moisture/MoistureMax;
            Mg=     W/obj.sub.VWCsat;  
                       
            % Eq. 16 - Sailor 2008... the second part is from EPlus code
            Lef = 1.91846e6*(Tf/(Tf-33.91))^2;
            if(obj.T_plants < 273.15) % Less than 0 C
                Lef = 2.838e6;
            end
            
            % Saturation pressure.. Eq. 29 and 30 -- Sailor.
            esf = 611.2*exp(17.67*(Tf-273.15)/(Tf-29.65));            
            qsf = 0.622*esf/(Pa-esf); % DIFFERENCE WITH E+... bug?
            esg = 611.2*exp(17.67*(Tg-273.15)/(Tg-29.65)); %Pa saturation vapor pressure
			qsg = 0.622*esg/(Pa-esg); %Saturation mixing ratio at ground surface temperature.
            
            
            % Eq. 15 - Sailor 2008
            qaf =((1-sigmaf)*qa+sigmaf*(0.3*qa+0.6*qsf*rn+0.1*qsg*Mg))/(1-sigmaf*(0.6*(1.0-rn)+0.1*(1.0-Mg)));
			
            
            % Eq. 14 - Sailor 2008            
            Lf = Lef*obj.plant.LAI*Rhoaf*Cf*Waf*rn*(qaf-qsf);
			
            if(R_sh > 200)
                1+1;
            end
            
            %% Soil energy budget            
            %%%%%%%%%%%%%%%%%%%%%
             
                                    
            % Derivative of Saturation vapor pressure, which is used in the
            % calculation of derivative saturation specific humidity
            %Desf = 611.2*exp(17.67*((Tf-273.15)/(Tf-29.65)))*(17.67*(Tf-273.15)*(-1.0)*(Tf- 29.65)^(-2)+17.67/(Tf-29.65));
			%dqf = ( ( 0.622 * Pa ) / ( Pa - esf )^2 ) * Desf; %Derivative of saturation specific humidity
			
            %Latent heat vaporiobj.plant.ztion  at the ground temperature
			Leg = 1.91846e6*(Tg/(Tg-33.91))^2;
			%Check to see if ice is sublimating or frost is forming.
            if ( obj.T_substrate < 273.15 ) 
                Leg = 2.838e6; %per FASST documentation p.15 after eqn. 37.
            end
            
            %Desg = 611.2*exp(17.67*((Tg-273.15)/(Tg-29.65)))*(17.67*(Tg-273.15)*(-1.0)*(Tg- 29.65)^(-2) + 17.67/(Tg-29.65));
			%dqg = (0.622*Pa/(Pa-esg)^2)*Desg;
            
            %Final Ground Atmosphere Energy Balance
			%Density of air at the soil surface temperature
			%Rhog = Pa/(air.R*Tg);
            Rhog = air.density(Pa,Tg);
            
            % Eq. 19 - Sailor 2008
            Rhoag=(Rhoa+Rhog)/2;
            
            % Eq. 23 - Sailor 2008
			Rib=2*consts.g*obj.plant.z*(Taf-Tg)/((Taf+Tg)*Waf*Waf); %Richardson Number
            
			% Eq. 22 - Sailor 2008... as in the EPlus code (they are not
			% the same)
            if ( Rib < 0.0 ) 
				Gammah = ( 1.0 - 16.0 * Rib)^(-0.5);
            else 
                if ( Rib >= 0.19 ) 
					Rib = 0.19;
                end
				Gammah = ( 1.0 - 5.0 * Rib)^(-0.5 ); % This 0.5 is added by EPlus... but is absent in Sailor's paper
            end
            
            % Eq. 21 - Sailor 2008
            Chng = (air.Kv/log(obj.plant.z/obj.plant.Zog))^2 / rch; % bulk transfer coefficient near ground
			
            % Eq. 20 - Sailor 2008
            Chg = Gammah*((1-sigmaf)*Chng+sigmaf*Chnf);
			
            % Eq. 18 - Sailor 2008
            sheatg = e0 + Rhoag * air.Cp * Chg * Waf; % added the e0 windless correction
			Hg = sheatg * ( Taf - Tg ); % sensible flux TO soil (W/m^2) DJS Jan 2011 (eqn. 32 in Frankenstein 2004)

            % Eq. 20 - Sailor 2008
			Chne = ( air.Kv / log( obj.plant.z / obj.plant.Zog ) )^2 / rche;
			
            % Latent heat flux in soil layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Eq. 26 - Sailor 2008
            Ce = Gammah*((1-sigmaf)*Chne + sigmaf * Chnf); % this is in fact Ceg in eq (28)

			% Eq. 25 - Sailor 2008
            qg = Mg * qsg + ( 1.0 - Mg ) * qaf; 
            
            % Eq. 24 - Sailor 2008... with an extra Mg?
            Lg = Ce * Leg * Waf * Rhoag * ( qaf - qg ) * Mg; 
            % In the FASST documentation there is NO Mg. However, in looking
			           
            
            % Conduction
            r_roof = obj.roof.depth/obj.roof.k;
            r_substrate = obj.sub.depth/obj.sub.thermalConductivity(W);
            Qcond=(TT(2)-obj.T_interior)/(r_roof + r_substrate + consts.rsi_roof);
            
            
            % Eq. 1 - Sailor
            Ff = sigmaf*(...
                R_sh*(1-obj.plant.rho) ... % short wave radiation
              + obj.plant.em*Latm ... % Long wave radiation from sky
              - obj.plant.em*consts.SB*Tf^4 ... % IR radiation loss
              ) ...
              + sigmaf*obj.sub.em*obj.plant.em*consts.SB*(Tg^4-Tf^4)/e_1 ... % IR radiation between plant and substrate
              + Hf ... % Sensible convective heat
              + Lf; % Latent heat (i.e. transpiration)
             
          % Eq. 17 - Sailor 2008 as shown in Eq. 31 of main report
          Fg = (1-sigmaf)*( ...
                R_sh*(1-obj.sub.rho)... % Short wave radiation
              + obj.sub.em*Latm ... % Ir radiation from sky
              - obj.sub.em*consts.SB*Tg^4 ...  % Emission from ground to sky
              )...              
              - sigmaf*obj.sub.em*obj.plant.em*consts.SB*(Tg^4-Tf^4)/e_1 ... % IR radiation between plant and substrate
              + Hg... % Sensible heat (i.e. convection)
              + Lg ... % Latent heat (i.e. evaporation)
              - Qcond; % Heat loss by conduction to ground
                  
            Res=[Ff Fg] ;           
            
        end % end resFUN
        
        function obj = moveForward(obj,data)            
            options=optimset('Display','off','TolFun',1e-6);
            
            Tguess = [obj.T_plants, obj.T_substrate];   
            
            Tss = (fsolve(@(TT) obj.ResFUN(TT,data),Tguess,options));
            
            obj.T_plants=(Tss(1));
            obj.T_substrate=(Tss(2)); 
            obj = obj.update(data);            
        end
        
        function T = getTemperature(obj, depth)
            %Film resistance coefficient
            consts = Constants;
            rsi_roof = consts.rsi_roof;
            
            %Roof resistance
            r_roof = obj.roof.k/obj.roof.depth;
            
            %Substrate resistance
            r_substrate = obj.sub.thermalConductivity(obj.VWC)/obj.sub.depth;
            
            %Total resistance
            total_r = r_roof + r_substrate + rsi_roof;
            
            %Total temperature difference
            total_dT = obj.T_interior - obj.T_substrate ;
            aux = total_dT/total_r;
            
            dT_roof = r_roof*aux;
            dT_sub = r_substrate*aux;
            
            % in this case, just interpolation
           if depth <= obj.sub.depth
               % in the substrate... interpolate
               T = obj.T_substrate + (depth/obj.sub.depth)*dT_sub;
           elseif depth <= (obj.sub.depth + obj.roof.depth)
               % in the roof
               depth = depth < obj.sub.depth; %translate to "roof" coordinates                
               T = obj.T_substrate + dT_sub + (depth/obj.roof.depth)*dT_roof;
           else
               % inside the space
               T = obj.T_interior;
           end
        end
        
        
        
    end % end methods section
    
end % end class