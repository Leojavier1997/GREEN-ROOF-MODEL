classdef SailorThermalMass
    properties
                
         %% Current state
        T_plants
        T_substrate
        T_interior
        VWC % volumetric water content
        dt = 60; %seconds between measurements
        interface_heat_flux % Heat flux at the interface                
        innerT
        
        % Plants balance
        plant_absorbed_solar
        plant_absorbed_ir_sky        
        transpiration
        plant_convection   
        
        % Both
        Qir                
                
        % Substrate balance
        substrate_solar_radiation
        substrate_infrared_radiation
        substrate_convection
        evaporation
        substrate_conduction 
        
        et_mm_hour  
        heating_load
        rs
        fsolar
        fvpd
        fvwc
        ftemp
        ra
        
        %% Materials
        roof
        plant
        sub
                
        
        %% Matrices for thermal mass
        C
        K
        
        %% other                
        L1
        n_roof_nodes % nnodes on the roof and other nnodes in the substrate
        n_sub_nodes
        samples
    end
    
    methods
        function obj = SailorThermalMass(plant,sub,roof,n_substrate_layers,n_roof_layers,Area,dt)
            
            obj.dt = dt*60;
            
            obj.roof = roof;
            obj.plant = plant;
            obj.sub = sub;
            
            %% Initialize
            obj.T_plants = 285; % K
            obj.T_substrate = 291; % K
            obj.T_interior = 293;
            
            
            %% Other      
            obj.n_sub_nodes = n_substrate_layers;
            obj.n_roof_nodes = n_roof_layers;
            obj.L1 = sqrt(Area);            
                 
            dxSubstrate = sub.depth/obj.n_sub_nodes;
            dxSupport = roof.depth/obj.n_roof_nodes;
            
            
            obj.samples = [0 dxSubstrate/2:dxSubstrate:(sub.depth-dxSubstrate/2) (sub.depth+dxSupport/2):dxSupport:(roof.depth+sub.depth-dxSupport/2) sub.depth+roof.depth];            
            
            %% Initialize inner temperatures by interpolating
            obj.innerT = interp1([0 (sub.depth+roof.depth)],[obj.T_substrate obj.T_interior],obj.samples(2:end-1))';           
            
            %% Create matrices
            %obj = obj.setMatrices();
            
        end
        
        function obj = setMatrices(obj)
            %% Update the matrices
            
             %Create empty matrix
            total_nodes = obj.n_roof_nodes + obj.n_sub_nodes;
            obj.C = zeros(total_nodes,total_nodes);
            obj.K = zeros(total_nodes,total_nodes);
            
            % Define parameters            
            n = 1; %Helper
            
            % Roof properties do not change
            mcRoof = (obj.roof.density*obj.roof.depth/obj.n_roof_nodes)*obj.roof.Cp;            
            rRoof = obj.roof.depth/obj.roof.k/obj.n_roof_nodes;
            
            % Connect before substrate (connection with T_substrate)
            rSub = obj.sub.depth/(obj.sub.thermalConductivity(obj.VWC(1))*obj.n_sub_nodes);
            obj.K(1,1) = 2/rSub;
            
            % Connect within the substate   
            for i=1:obj.n_sub_nodes-1
                mcSub = (obj.sub.depth/obj.n_sub_nodes)*obj.sub.rhoCp(obj.VWC(n));          
                rSub1 = obj.sub.depth/(obj.sub.thermalConductivity(obj.VWC(n))*obj.n_sub_nodes);
                rSub2 = obj.sub.depth/(obj.sub.thermalConductivity(obj.VWC(n+1))*obj.n_sub_nodes);
            
                obj.C(n,n) = mcSub;
                obj.K(n:n+1,n:n+1)=obj.K(n:n+1,n:n+1)+[1,-1;-1,1]./(rSub1/2+rSub2/2);
                n = n+1;
            end
            
            % Connect interface between both
            rSub = obj.sub.depth/(obj.sub.thermalConductivity(obj.VWC(n))*obj.n_sub_nodes);
            obj.C(n,n) = mcSub;
            obj.K(n:n+1,n:n+1)=obj.K(n:n+1,n:n+1)+[1,-1;-1,1]./(rSub/2 + rRoof/2);
            n=n+1;
            
            % Connect within the roof            
            for i=1:obj.n_roof_nodes-1
                obj.C(n,n) = mcRoof;
                obj.K(n:n+1,n:n+1)=obj.K(n:n+1,n:n+1)+[1,-1;-1,1]./rRoof;
                n = n+1;
            end
            
            % Connect final one.
            obj.C(n,n) = mcRoof;
            consts = Constants;
            obj.K(n,n) = obj.K(n,n)+1/(consts.rsi_roof + rRoof/2);
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % *MODIFICACIÓN 1 VELOCIDAD MÁXIMA DEL AIRE*
            %%%%%original%%%%%
            U=      max(3,data.U);       % Wind speed [m/s]
            %%%%%PRUEBA 1 SIN MÁX%%%%%
            %U=max(0.001,data.U); %wind speed 
            %U=max(2,data.U);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            W=      obj.VWC(1);    %VWC            
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
            
            display("Taf "+Taf+"Tair "+Tair+"Tf "+ Tf +"Tg " +Tg)
            
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
            obj.ra = ra;
            
            % Eq. 11a - Sailor 2008
            f1 = 1/min(1,(0.004*R_sh+0.005)/(0.81*(0.004*R_sh+1.0)));
            
            % Unsure about these definitions... made to fit EPlus with
            MeanRootMoisture = obj.getVWC(obj.sub.depth/2);
            MoistureResidual = obj.sub.VWCresidual;
            MoistureMax = obj.sub.VWCsat;
            %display(MoistureResidual + " " +MoistureMax)
            
            % Eq. 11b - Sailor 2008... f2 does not become zero. This idea
            % comes from EnergyPlus code            
            if(MeanRootMoisture <= MoistureResidual)
                f2 = 1000;
            else
                f2 = (MoistureMax-MoistureResidual)/(MeanRootMoisture-MoistureResidual);                
            end
            
            % Eq. 11c - Sailor 2008
            f3 = 1; %1/exp(-0*(esf-eair));
            %Ecuación original Sailor:
            %r_s = obj.plant.rsmin*f1*f2*f3/obj.plant.LAI;            
            
            e_s  = 610.8*exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000;     
            e_sf = 610.8*exp(17.27*(TT(1)-273.15)/(TT(1)-273.15+237.3))/1000;
            e_air= e_s*RH;
            VPD_f=e_sf-e_air;
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Ecuación Resistencia estomática RLM:

            %Modelo 5 or New Linear Model: rs_min, Rsh, VWC, LAI, VPD
            r_s=184.18+1.987*obj.plant.rsmin-0.157*R_sh-247.116*MeanRootMoisture-34.592*obj.plant.LAI+21.568*VPD_f;
            if ( r_s <obj.plant.rsmin) 
                r_s=   obj.plant.rsmin;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         

            obj.rs = r_s;
            
            obj.fsolar=f1;
            obj.fvpd=f3;
            obj.fvwc=f2;
            obj.ftemp=1;
            
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
            
            
            
            %% Update inner temperatures
            % find temperatures within the subtrate/concrete
            rSub = (obj.sub.depth/obj.sub.thermalConductivity(obj.VWC(1)))/obj.n_sub_nodes;
            rRoof = (obj.roof.depth/obj.roof.k)/obj.n_roof_nodes;
            
            f      = zeros(length(obj.innerT),1);
            f(1)   = 2*obj.T_substrate/rSub;
            f(end) = obj.T_interior/(consts.rsi_roof+rRoof/2);
            %obj.innerT = inv(obj.C/obj.dt + obj.K)*(obj.C/obj.dt*obj.innerT + f);
            %obj.innerT = obj.innerT + obj.C\(f - obj.K*obj.innerT)*obj.dt;
            total_nodes = obj.n_roof_nodes + obj.n_sub_nodes;
            R = eye(total_nodes,total_nodes)+(obj.C\(obj.dt*obj.K)); 
            obj.innerT = R\(obj.innerT+(obj.C\f*obj.dt));  
            Qcond = 2*(TT(2) - obj.innerT(1))/rSub;
                                               
            %% UPDATE                 
            % Plants balance
            obj.plant_absorbed_solar = R_sh*(1-obj.plant.rho)*sigmaf;
            obj.plant_absorbed_ir_sky = (obj.plant.em*Latm - obj.plant.em*consts.SB*Tf^4)*sigmaf;
            obj.transpiration = Lf;
            obj.plant_convection = Hf;

            % Both
            obj.Qir = -sigmaf*obj.sub.em*obj.plant.em*consts.SB*(Tg^4-Tf^4)/e_1;
              
            % Substrate balance
            obj.substrate_solar_radiation = (1-sigmaf)*R_sh*(1-obj.sub.rho);
            obj.substrate_infrared_radiation = (1-sigmaf)*(obj.sub.em*Latm - obj.sub.em*consts.SB*Tg^4);
            obj.substrate_convection = Hg;
            obj.evaporation = Lg;
            obj.substrate_conduction = -Qcond;                                                 
            
            
            obj.et_mm_hour = -(Lg/Leg + Lf/Lef)*3600; 
            
            
            %% Update interface_heat_flux
            rSub = (obj.sub.depth/obj.sub.thermalConductivity(obj.VWC(obj.n_sub_nodes)))/obj.n_sub_nodes; % Last node            
            deltaT = obj.innerT(obj.n_sub_nodes)-obj.innerT(obj.n_sub_nodes+1); %one node on each.  
            Qcond = 2*deltaT/(rSub+rRoof);
            obj.interface_heat_flux = Qcond;
            
            Tlosa = obj.getTemperature(obj.sub.depth + obj.roof.depth);
            obj.heating_load = ( obj.T_interior - Tlosa)/consts.rsi_roof;
             
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % *MODIFICACIÓN 1 VELOCIDAD MÁXIMA DEL AIRE*
            %%%%%original%%%%%
            U=      max(3,data.U);       % Wind speed [m/s]
            %%%%%PRUEBA 1 SIN MÁX%%%%%
            %U=max(0.001,data.U); %wind speed 
            %U=max(2,data.U);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            W=      obj.VWC(1);    %VWC            
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
            MeanRootMoisture = obj.getVWC(obj.sub.depth/2);
            MoistureResidual = obj.sub.VWCresidual;
            MoistureMax = obj.sub.VWCsat;
            
            % Eq. 11b - Sailor 2008... f2 does not become zero. This idea
            % comes from EnergyPlus code            
            if(MeanRootMoisture <= MoistureResidual)
                f2 = 1000;
            else
                f2 = (MoistureMax-MoistureResidual)/(MeanRootMoisture-MoistureResidual);                
            end
            
            
            % Eq. 11c - Sailor 2008
            f3 = 1; %1/exp(0*(esf-air));
            
            %Ecuación original Sailor:
            %r_s = obj.plant.rsmin*f1*f2*f3/obj.plant.LAI;            

            e_s  = 610.8*exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000;     
            e_sf = 610.8*exp(17.27*(TT(1)-273.15)/(TT(1)-273.15+237.3))/1000;
            e_air= e_s*RH;
            VPD_f=e_sf-e_air;
            


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Ecuación Resistencia estomática RLM:
   
            %Modelo 5 or New Linear Model: rs_min, Rsh, VWC, LAI, VPD
            r_s=184.18+1.987*obj.plant.rsmin-0.157*R_sh-247.116*MeanRootMoisture-34.592*obj.plant.LAI+21.568*VPD_f;
            if ( r_s <obj.plant.rsmin) 
                r_s=   obj.plant.rsmin;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            %display("R_sh " + R_sh +" MeanRootMoisture " + MeanRootMoisture + " Rs_min " + obj.plant.rsmin +" Tair" + Tair +" rs "+r_s);
            
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
			           
            
            %% CONDUCTION
                                    
            % eq. 8 Tabares modified to only consider the conduction
            % through the first layer of substrate. 
            
            % find temperatures within the subtrate/concrete            
            rSub = obj.sub.depth/obj.sub.thermalConductivity(obj.VWC(1))/obj.n_sub_nodes;
            rRoof = obj.roof.depth/obj.roof.k/obj.n_roof_nodes;
            
            f = zeros(length(obj.innerT),1);
            f(1) =  2*TT(2)/rSub;
            f(end)= obj.T_interior/(consts.rsi_roof+rRoof/2);
           
            %iT = inv(obj.C/obj.dt - obj.K)*((obj.C/obj.dt)*obj.innerT - f);
            %iT = obj.innerT + obj.C\(f - obj.K*obj.innerT)*obj.dt;
            total_nodes = obj.n_roof_nodes + obj.n_sub_nodes;
            R = eye(total_nodes,total_nodes)+(obj.C\(obj.dt*obj.K)); 
            iT = R\(obj.innerT+(obj.C\f*obj.dt));  
            Qcond = 2*(TT(2) - iT(1))/rSub;
            
            
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
            
            obj = obj.setMatrices();
            
            options=optimset('Display','off','TolFun',1e-6);
            Tguess = [obj.T_plants, obj.T_substrate];              
            Tss = (fsolve(@(TT) obj.ResFUN(TT,data),Tguess,options));
            
            obj.T_plants=(Tss(1));
            obj.T_substrate=(Tss(2)); 
            obj = obj.update(data);            
        end
        
        function T = getTemperature(obj, depth)
            if depth > (obj.roof.depth + obj.sub.depth)                
                T = obj.T_interior;
            else
                consts = Constants;
                Rroof = obj.roof.depth/obj.roof.k/obj.n_roof_nodes/2;
                Ts=(Rroof*obj.T_interior+consts.rsi_roof*obj.innerT(end))/(consts.rsi_roof+Rroof);
                T = interp1([obj.samples],[obj.T_substrate obj.innerT' Ts],depth);                
            end
        end
        
        
        function vwc = getVWC(obj, depth)
            dx = (obj.sub.depth)/obj.n_sub_nodes;
                        
            if depth > obj.sub.depth                
                vwc = 0;
            else
                vwc = interp1([dx/2:dx:(obj.sub.depth-dx/2)],[obj.VWC],depth);                
            end
        end
        
        
        
    end % end methods section
    
end % end class