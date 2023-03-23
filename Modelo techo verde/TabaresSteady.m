classdef TabaresSteady
    properties
        %% Current state
        T_plants
        T_substrate
        T_interior        
        VWC % volumetric water content       
        interface_heat_flux
        evaporation
        transpiration
        substrate_convection
        plant_convection
        
        et_mm_hour        
        
        dt
        
        %% Materials
        roof
        plant
        sub
        
        %% other                
        L1
        samples % points where temperature is known.
    end
    
    methods
        function obj = TabaresSteady(plant,sub,roof,n_layers,Area,dt)
            
            obj.roof = roof;
            obj.plant = plant;
            obj.sub = sub;
            
            %% Initialize
            obj.T_plants = 285; % K
            obj.T_substrate = 291; % K
            obj.T_interior = 293;
            
            %% Other                        
            obj.L1 = sqrt(Area);            
            obj.VWC = sub.VWC;
            obj.samples = [0,sub.depth,sub.depth+roof.depth];
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
            
            T=[obj.T_plants,obj.T_substrate];
            
            %% Derived inputs
            Mu=air.kinematicViscosity(Tair);   % ***** viscosidad cinematica
            dens_a= air.density(Pa,Tair);   % density air            
            Tfilm=  (Tair+T(1))/2;
            Beta=   1/Tfilm;
            
            
            %% OTHER            
            % From table 5.
            kpor  = obj.sub.phi*air.k+(1-obj.sub.phi)*obj.plant.k;            
            % Auxiliar variable
            Mg=     W/obj.sub.VWCsat;        
            % No idea where this comes from.
            alphapor= kpor/(dens_a*air.Cp);
            
            % eq. 6 Tabares ... constants may vary?
            rsub=   34.52*Mg^(-3.2678);
            
            
            %%   DIMENSIONLESS NUMBERS
            % Reynolds **note* with L1 not L
            Re=     dens_a*U*obj.L1/Mu;
            % Prandtl
            Pr=     air.Cp*Mu/air.k;
            % Grashof
            Gr=     abs(consts.g*Beta*dens_a^2*(T(1)-Tair)*(obj.L1^3)/(Mu^2));
            % Raleigh
            Ra=     Gr*Pr;
            % Lewis
            Le=     1;
            
            % Nusselt. eq. 4 
            if (Gr<(0.068*Re^2.2))
                Nu= 3+1.25*0.025*Re^0.8;
            elseif (Gr>(0.068*Re^2.2))&&(Gr<(55.3*Re^(5/3)))
                Nu= 2.7*((Gr/(Re^2.2))^(1/3))*(3*(15/4)+(15/16)*0.0253*Re^0.8);
            else
                Nu= 0.15*(Ra^(1/3));
            end
            
            % available in Table 5. Tabares
            Pe=     0.3*obj.L1*U/alphapor; % CHECK
            hpor  = kpor*1.128*Pe^0.5/obj.L1; % CHECK            
            hconv= 1.5*Nu*air.k/obj.L1; % plantBB was 1.5 then 15
            %ras =  dens_a*air.Cp*(Le^(2/3))/hconv;                   
            ra  =  dens_a*air.Cp*(Le^(2/3))/hconv;  %ORIGINAL
            hsub =  hpor*hconv/(hpor+hconv); % Defined in texti, below Eq. 19
            
            % Vapor Pressures [kPa]
            
            %e_s  = 610.8*exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000;     
            e_s = 611.2*exp(17.67*(Tair-273.15)/(Tair-29.65))/1000;
            %e_sf = 610.8*exp(17.27*(T(1)-273.15)/(T(1)-273.15+237.3))/1000;
            e_sf = 611.2*exp(17.67*(T(1)-273.15)/(T(1)-29.65))/1000;
            %e_ss = 610.8*exp(17.27*(T(2)-273.15)/(T(2)-273.15+237.3))/1000;   
            e_ss = 611.2*exp(17.67*(T(2)-273.15)/(T(2)-29.65))/1000;
            e_air= e_s*RH;

                                                          
            %% SHORT WAVE RADIATION
            
            % eq. 12 Tabares
            Rsh_f = (1-obj.plant.foliage_rho-obj.plant.tau_fsol)*(1+obj.plant.tau_fsol*obj.sub.rho)*R_sh;
            
            % eq. 13 Tabares
            Rsh_s = obj.plant.tau_fsol*(1-obj.sub.rho)*R_sh;
            
            %% LONG WAVE RADIATION
            
            % eq. 14 Tabares
            Qir_f= (1-obj.plant.tau_fir)*obj.plant.em*consts.SB*(T(1)^4-Tsky^4);
            
            % eq. 15 Tabares
            Qir_scov=   (obj.plant.tau_fir)*obj.sub.em*consts.SB*((T(2)^4-Tsky^4));
            
            % eq. 17 Tabares
            em_1 = (1/obj.sub.em)+(1/obj.plant.em)-1;
            Qir_sp = (1-obj.plant.tau_fir)*consts.SB*(T(1)^4-T(2)^4)/em_1;
            
            %% CONVECTION
            
            % eq. 18 Tabares -- 1.5*LAI*hconv*(Tplant-Tair) ... 1.5?
            obj.substrate_convection=   1.5* obj.plant.LAI*hconv*(T(1)-Tair);
            
            % eq. 19 Tabares
            obj.plant_convection=    hsub*(T(2)-Tair);
            
            %% EVAPOTRANSPIRATION
            
            % eq. 22 Tabares
            f_sol=  1+exp(-0.034*(R_sh-3.5));
            
            % eq. 24 Tabares... extended?
            %f_VPD=  1/(1-0.41*log(e_sf-e_air));
            VPD=e_ss-e_air;
            VPD_f=e_sf-e_air;

            if VPD_f > 0
                f_VPD=  (1-0.41*log(VPD_f));
            else
                f_VPD= 1;
            end
            if f_VPD>1
                f_VPD=1;
            elseif f_VPD<0
                f_VPD=0.05;
            end

            % eq. 25 Tabares
            f_temp= 1/(1-0.0016*(35-(T(1)-273.15))^2);
            
            % eq. 23 Tabares
            if (W>0.7*obj.sub.VWCfc) % W_fc vs W_sat
                f_W=    1;
            elseif (W<=0.7*obj.sub.VWCfc) && (W>obj.sub.VWCresidual)
                f_W=    (0.7*obj.sub.VWCfc-obj.sub.VWCresidual)/(W-obj.sub.VWCresidual);
            else
                f_W=    1000;
            end
            
            % eq. 21 Tabares
            f_hum=1/(f_VPD);    
            rs = obj.plant.rsmin*f_sol*f_hum*f_W*f_temp/obj.plant.LAI;

            % eq. 20 Tabares
            %Qt= obj.plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ras)); 
            obj.transpiration= obj.plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ra)); 
            
            % eq. 5 Tabares
            %evaporation= dens_a*air.Cp*(e_s-e_air)/(air.gamma*(rsub+ras));    
            obj.evaporation= dens_a*air.Cp*VPD/(air.gamma*(rsub+ra));    
               
            % Eq. 16 - Sailor 2008... the second part is from EPlus code
            Tg = obj.T_substrate;
            Tf = obj.T_plants;
            Lef = 1.91846e6*(Tf/(Tf-33.91))^2;
            if(obj.T_plants < 273.15) % Less than 0 C
                Lef = 2.838e6;
            end
            
            Leg = 1.91846e6*(Tg/(Tg-33.91))^2;
            if(obj.T_plants < 273.15) % Less than 0 C
                Lef = 2.838e6;
            end
            obj.et_mm_hour = -(obj.evaporation/Leg + obj.transpiration/Lef)*3600;   
            
            %% Update interface_heat_flux
            r_roof = obj.roof.depth/obj.roof.k;
            r_substrate = obj.sub.depth/obj.sub.thermalConductivity(W);            
            obj.interface_heat_flux = (obj.T_substrate-obj.T_interior)/(r_roof + r_substrate + consts.rsi_roof);
            
            if R_sh > 200
                1+1;
            end
            
        end
        
        function [Res] = ResFUN(obj,T,data)
            
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
            
            %% Derived inputs
            Mu=air.kinematicViscosity(Tair);   % ***** viscosidad cinematica
            dens_a= air.density(Pa,Tair);   % density air            
            Tfilm=  (Tair+T(1))/2;
            Beta=   1/Tfilm;
            
            
            %% OTHER            
            % From table 5.
            kpor  = obj.sub.phi*air.k+(1-obj.sub.phi)*obj.plant.k;            
            % Auxiliar variable
            Mg=     W/obj.sub.VWCsat;        
            % No idea where this comes from.
            alphapor= kpor/(dens_a*air.Cp);
            
            % eq. 6 Tabares ... constants may vary?
            rsub=   34.52*Mg^(-3.2678);
            
            
            %%   DIMENSIONLESS NUMBERS
            % Reynolds **note* with L1 not L
            Re=     dens_a*U*obj.L1/Mu;
            % Prandtl
            Pr=     air.Cp*Mu/air.k;
            % Grashof
            Gr=     abs(consts.g*Beta*dens_a^2*(T(1)-Tair)*(obj.L1^3)/(Mu^2));
            % Raleigh
            Ra=     Gr*Pr;
            % Lewis
            Le=     1;
            
            % Nusselt. eq. 4 
            if (Gr<(0.068*Re^2.2))
                Nu= 3+1.25*0.025*Re^0.8;
            elseif (Gr>(0.068*Re^2.2))&&(Gr<(55.3*Re^(5/3)))
                Nu= 2.7*((Gr/(Re^2.2))^(1/3))*(3*(15/4)+(15/16)*0.0253*Re^0.8);
            else
                Nu= 0.15*(Ra^(1/3));
            end
            
            % available in Table 5. Tabares
            Pe=     0.3*obj.L1*U/alphapor; % CHECK
            hpor  = kpor*1.128*Pe^0.5/obj.L1; % CHECK            
            hconv= 15*Nu*air.k/obj.L1; % plantBB was 1.5 then 15
            ra=dens_a*air.Cp*(Le^(2/3))/hconv;  %ORIGINAL
            hsub =  hpor*hconv/(hpor+hconv); % CHECK
            
            % Vapor Pressures [kPa]
            
             %e_s  = 610.8*exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000;     
            e_s = 611.2*exp(17.67*(Tair-273.15)/(Tair-29.65))/1000;
            %e_sf = 610.8*exp(17.27*(T(1)-273.15)/(T(1)-273.15+237.3))/1000;
            e_sf = 611.2*exp(17.67*(T(1)-273.15)/(T(1)-29.65))/1000;
            %e_ss = 610.8*exp(17.27*(T(2)-273.15)/(T(2)-273.15+237.3))/1000;   
            e_ss = 611.2*exp(17.67*(T(2)-273.15)/(T(2)-29.65))/1000;
            e_air= e_s*RH;

                                                          
            %% SHORT WAVE RADIATION
            
            % eq. 12 Tabares
            Rsh_f = (1-obj.plant.foliage_rho-obj.plant.tau_fsol)*(1+obj.plant.tau_fsol*obj.sub.rho)*R_sh;
            
            % eq. 13 Tabares
            Rsh_s = obj.plant.tau_fsol*(1-obj.sub.rho)*R_sh;
            
            %% LONG WAVE RADIATION
            
            % eq. 14 Tabares
            Qir_f= (1-obj.plant.tau_fir)*obj.plant.em*consts.SB*(T(1)^4-Tsky^4);
            
            % eq. 15 Tabares
            Qir_scov=   (obj.plant.tau_fir)*obj.sub.em*consts.SB*((T(2)^4-Tsky^4));
            
            % eq. 17 Tabares
            em_1 = (1/obj.sub.em)+(1/obj.plant.em)-1;
            Qir_sp = (1-obj.plant.tau_fir)*consts.SB*(T(1)^4-T(2)^4)/em_1;
            
            %% CONVECTION
            
            % eq. 18 Tabares -- 1.5*LAI*hconv*(Tplant-Tair) ... 1.5?
            substrate_convection=   1.5* obj.plant.LAI*hconv*(T(1)-Tair);
            
            % eq. 19 Tabares
            plant_convection=    hsub*(T(2)-Tair);
            
            %% EVAPOTRANSPIRATION
            
            % eq. 22 Tabares
            f_sol=  1+exp(-0.034*(R_sh-3.5));
            
            % eq. 24 Tabares... extended?
            %f_VPD=  1/(1-0.41*log(e_sf-e_air));
            VPD=e_ss-e_air;
            VPD_f=e_sf-e_air;

            if VPD_f > 0
                f_VPD=  (1-0.41*log(VPD_f));
            else
                f_VPD= 1;
            end
            if f_VPD>1
                f_VPD=1;
            elseif f_VPD<0
                f_VPD=0.05;
            end

            % eq. 25 Tabares
            f_temp= 1/(1-0.0016*(35-(T(1)-273.15))^2);
            
            % eq. 23 Tabares
            if (W>0.7*obj.sub.VWCfc) % W_fc vs W_sat
                f_W=    1;
            elseif (W<0.7*obj.sub.VWCfc) && (W>obj.sub.VWCresidual)
                f_W=    (0.7*obj.sub.VWCfc-obj.sub.VWCresidual)/(W-obj.sub.VWCresidual);
            else
                f_W=    1000;
            end
            
            % eq. 21 Tabares
            f_hum=1/(f_VPD);    
            rs = obj.plant.rsmin*f_sol*f_hum*f_W*f_temp/obj.plant.LAI;

            % eq. 20 Tabares
            Qt= obj.plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ra)); 
            
            % eq. 5 Tabares
            %evaporation= dens_a*air.Cp*(e_s-e_air)/(air.gamma*(rsub+ras));    
            evaporation= dens_a*air.Cp*VPD/(air.gamma*(rsub+ra));    
            
            %% CONDUCTION
            
            % eq. 8 Tabares modified to include the concrete below
            r_roof = obj.roof.depth/obj.roof.k;
            r_substrate = obj.sub.depth/obj.sub.thermalConductivity(W);
            Qcond=(T(2)-obj.T_interior)/(r_roof + r_substrate + consts.rsi_roof);
            
            %% ENERGY (& MASS) BALANCE
            
            % Foliage energy balance
            Ef= - Rsh_f + substrate_convection + Qir_f + Qt + Qir_sp; % Eq 14 Camilo
            
            % Substrate energy balance
            Es= - Rsh_s - Qir_sp + plant_convection + Qir_scov + evaporation + Qcond; % Eq 15 Camilo
                        
            %% RETURN
            Res=[Ef Es] ;
            
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