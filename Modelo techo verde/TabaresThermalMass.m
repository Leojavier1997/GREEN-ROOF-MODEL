classdef TabaresThermalMass
    properties
        %% Current state        
        T_plants
        T_substrate
        T_interior        
        VWC % volumetric water content     
        
        interface_heat_flux
        
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
        
        
        innerT
        dt
        
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
        function obj = TabaresThermalMass(plant,sub,roof,n_substrate_layers,n_roof_layers,Area,dt)
            
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
            obj.innerT = real(obj.innerT);
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
            %U=data.U; %wind speed 
            %U=max(2,data.U);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %W=      obj.VWC(1);    %VWC            
            
            T = [obj.T_plants,obj.T_substrate];
            
            %% Derived inputs
            Mu=air.kinematicViscosity(Tair);   % ***** viscosidad cinematica
            dens_a= air.density(Pa,Tair);   % density air            
            Tfilm=  (Tair+T(1))/2;
            Beta=   1/Tfilm;
            
            
            %% OTHER            
            % From table 5.
            kpor  = obj.sub.phi*air.k+(1-obj.sub.phi)*obj.plant.k;            
            % Auxiliar variable
            Mg=     obj.VWC(1)/obj.sub.VWCsat;        
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
            ras=  dens_a*air.Cp*(Le^(2/3))/hconv;                   
            ra=dens_a*air.Cp*(Le^(2/3))/hconv;  %ORIGINAL
            obj.ra = ra;
            hsub =  hpor*hconv/(hpor+hconv); % CHECK
            
            % Vapor Pressures [kPa]
            
            e_s  = 610.8*exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000;     
            e_sf = 610.8*exp(17.27*(T(1)-273.15)/(T(1)-273.15+237.3))/1000;
            e_ss = 610.8*exp(17.27*(T(2)-273.15)/(T(2)-273.15+237.3))/1000;   
            e_air= e_s*RH;

                                                          
            %% SHORT WAVE RADIATION
                
            % eq. 12 Tabares
            obj.plant_absorbed_solar = (1-obj.plant.foliage_rho-obj.plant.tau_fsol)*(1+obj.plant.tau_fsol*obj.sub.rho)*R_sh;
            
            % eq. 13 Tabares
            obj.substrate_solar_radiation = obj.plant.tau_fsol*(1-obj.sub.rho)*R_sh;
            
            %% LONG WAVE RADIATION
         
            % eq. 14 Tabares
            obj.plant_absorbed_ir_sky = (1-obj.plant.tau_fir)*obj.plant.em*consts.SB*(T(1)^4-Tsky^4)
            

            % eq. 15 Tabares
            obj.substrate_infrared_radiation =   -(obj.plant.tau_fir)*obj.sub.em*consts.SB*((T(2)^4-Tsky^4));
            
            % eq. 17 Tabares
            em_1 = (1/obj.sub.em)+(1/obj.plant.em)-1;
            obj.Qir = (1-obj.plant.tau_fir)*consts.SB*(T(1)^4-T(2)^4)/em_1;
            
            %% CONVECTION                            
            % eq. 18 Tabares -- 1.5*LAI*hconv*(Tplant-Tair) ... 1.5?
            obj.plant_convection=   -1.5* obj.plant.LAI*hconv*(T(1)-Tair);
            
            % eq. 19 Tabares
            obj.substrate_convection=    -hsub*(T(2)-Tair);
            
                
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
            if f_temp < 0
                f_temp = 10000;
            end
            % eq. 23 Tabares
            rootW = obj.getVWC(obj.sub.depth/2);
            if ( rootW >0.7*obj.sub.VWCfc) % W_fc vs W_sat
                f_W=    1;
            elseif (rootW <= 0.7*obj.sub.VWCfc) && (rootW > obj.sub.VWCresidual)
                f_W=    (0.7*obj.sub.VWCfc-obj.sub.VWCresidual)/(rootW-obj.sub.VWCresidual);
            else
                f_W=    1000;
            end
            
            %display(obj.sub.VWCresidual + " " +obj.sub.VWCfc)
            % eq. 21 Tabares
            %rs = obj.plant.rsmin*f_sol*f_VPD*f_W*f_temp/obj.plant.LAI;
            f_hum=1/(f_VPD);    
            %Ecuación original
            %rs = (obj.plant.rsmin/obj.plant.LAI)*f_sol*f_hum*f_W*f_temp;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Ecuación Resistencia estomática RLM:

            %Modelo 5 or New Linear Model: rs_min, Rsh, VWC, LAI, VPD
            rs=184.18+1.987*obj.plant.rsmin-0.157*R_sh-247.116*rootW-34.592*obj.plant.LAI+21.568*VPD_f;
            if ( rs <obj.plant.rsmin) 
                rs=   obj.plant.rsmin;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.rs = rs;
            obj.fsolar=f_sol;
            obj.fvpd=f_hum;
            obj.fvwc=f_W;
            obj.ftemp=f_temp;
            
            % eq. 20 Tabares
            %Qt= obj.plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ras)); 
            obj.transpiration= -obj.plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ra)); 
            
            % eq. 5 Tabares
            %evaporation= dens_a*air.Cp*(e_s-e_air)/(air.gamma*(rsub+ras));    
            obj.evaporation= -dens_a*air.Cp*VPD/(air.gamma*(rsub+ras));                               
        
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
            display("Tair "+Tair+"Tf "+ Tf +"Tg " +Tg)

            
                                             
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
            obj.innerT  = R\(obj.innerT+(obj.C\f*obj.dt)); 
            
            Qcond = 2*(T(2) - obj.innerT(1))/rSub;              
            obj.substrate_conduction = -Qcond;
            
            %% Update interface_heat_flux            
            rSub = (obj.sub.depth/obj.sub.thermalConductivity(obj.VWC(obj.n_sub_nodes)))/obj.n_sub_nodes; % Last node
            deltaT = obj.innerT(obj.n_sub_nodes)-obj.innerT(obj.n_sub_nodes+1); %one node on each.            
            obj.interface_heat_flux = 2*deltaT/(rSub+rRoof);
            
            Tlosa = obj.getTemperature(obj.sub.depth + obj.roof.depth);
            obj.heating_load = ( obj.T_interior - Tlosa)/consts.rsi_roof;
            
            
        end
        
        function [Res] = ResFUN(obj,T,data)
            obj.innerT = real(obj.innerT);
           
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
            %U=data.U; %wind speed 
            %U=max(2,data.U);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %W=      obj.VWC(1);    %VWC            
            
            %% Derived inputs
            Mu=air.kinematicViscosity(Tair);   % ***** viscosidad cinematica
            dens_a= air.density(Pa,Tair);   % density air            
            Tfilm=  (Tair+T(1))/2;
            Beta=   1/Tfilm;
            
            
            %% OTHER            
            % From table 5.
            kpor  = obj.sub.phi*air.k+(1-obj.sub.phi)*obj.plant.k;            
            % Auxiliar variable
            Mg=     obj.VWC(1)/obj.sub.VWCsat;        
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
            ras=  dens_a*air.Cp*(Le^(2/3))/hconv;                   
            ra=dens_a*air.Cp*(Le^(2/3))/hconv;  %ORIGINAL
            hsub =  hpor*hconv/(hpor+hconv); % CHECK
            
            % Vapor Pressures [kPa]
            
            e_s  = 610.8*exp(17.27*(Tair-273.15)/(Tair-273.15+237.3))/1000;     
            e_sf = 610.8*exp(17.27*(T(1)-273.15)/(T(1)-273.15+237.3))/1000;
            e_ss = 610.8*exp(17.27*(T(2)-273.15)/(T(2)-273.15+237.3))/1000;   
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
            plant_convection=   1.5* obj.plant.LAI*hconv*(T(1)-Tair);
            
            % eq. 19 Tabares
            substrate_convection=    hsub*(T(2)-Tair);
            
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
            if f_temp < 0
                f_temp = 10000;
            end
            
            % eq. 23 Tabares
            rootW = obj.getVWC(obj.sub.depth/2);
            if (rootW > 0.7*obj.sub.VWCfc) % W_fc vs W_sat
                f_W=    1;
            elseif (rootW < 0.7*obj.sub.VWCfc) && (rootW > obj.sub.VWCresidual)
                f_W=    (0.7*obj.sub.VWCfc-obj.sub.VWCresidual)/(rootW-obj.sub.VWCresidual);
            else
                f_W=    1000;
            end
            
            
            % eq. 21 Tabares
            %rs = obj.plant.rsmin*f_sol*f_VPD*f_W*f_temp/obj.plant.LAI;
            f_hum=1/(f_VPD);   
            %Ecuación original
            %rs = (obj.plant.rsmin/obj.plant.LAI)*f_sol*f_hum*f_W*f_temp;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Ecuación Resistencia estomática RLM:


            %Modelo 5 or New Linear Model: rs_min, Rsh, VWC, LAI, VPD
            rs=184.18+1.987*obj.plant.rsmin-0.157*R_sh-247.116*rootW-34.592*obj.plant.LAI+21.568*VPD_f;
            if ( rs <obj.plant.rsmin) 
                rs=   obj.plant.rsmin;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % eq. 20 Tabares
            %Qt= obj.plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ras)); 
            Qt= obj.plant.LAI*dens_a*air.Cp*(e_sf-e_air)/(air.gamma*(rs+ra)); 
            
            % eq. 5 Tabares
            %evaporation= dens_a*air.Cp*(e_s-e_air)/(air.gamma*(rsub+ras));    
            evaporation= dens_a*air.Cp*VPD/(air.gamma*(rsub+ras));    
            
            %% CONDUCTION
                                    
            % eq. 8 Tabares modified to only consider the conduction
            % through the first layer of substrate. 
            
            % find temperatures within the subtrate/concrete            
            rSub = obj.sub.depth/obj.sub.thermalConductivity(obj.VWC(1))/obj.n_sub_nodes; % Using thermal conductivity of the first node.
            rRoof = obj.roof.depth/obj.roof.k/obj.n_roof_nodes;
            
            f = zeros(length(obj.innerT),1);
            f(1) =  2*T(2)/rSub;
            f(end)= obj.T_interior/(consts.rsi_roof+rRoof/2);
           
            %iT = inv(obj.C/obj.dt - obj.K)*((obj.C/obj.dt)*obj.innerT - f);
            %iT = obj.innerT + obj.C\(f - obj.K*obj.innerT)*obj.dt;
            total_nodes = obj.n_roof_nodes + obj.n_sub_nodes;
            R = eye(total_nodes,total_nodes)+(obj.C\(obj.dt*obj.K)); 
            iT = R\(obj.innerT+(obj.C\f*obj.dt)); 
            
            Qcond = 2*(T(2) - iT(1))/rSub;
            
            %% ENERGY (& MASS) BALANCE
            
            % Foliage energy balance
            Ef= - Rsh_f + plant_convection + Qir_f + Qt + Qir_sp; % Eq 14 Camilo
            
            % Substrate energy balance
            Es= - Rsh_s - Qir_sp + substrate_convection + Qir_scov + evaporation + Qcond; % Eq 15 Camilo
                    
           
            %% RETURN
            Res=[Ef Es] ;
            
        end % end resFUN
        
        
        function obj = moveForward(obj,data)
            
            
            obj = obj.setMatrices();
            
            options=optimset('Display','off','TolFun',1e-6);
            Tguess = [obj.T_plants obj.T_substrate];                                                      
            Tss = (fsolve(@(TT) obj.ResFUN(TT,data),Tguess,options));
            obj.T_plants=(Tss(1));
            obj.T_substrate=(Tss(2)); 
            obj = obj.update(data);   
            
            obj.T_plants = real(obj.T_plants);
            obj.T_substrate = real(obj.T_substrate);
            obj.interface_heat_flux = real(obj.interface_heat_flux);
            obj.interface_heat_flux = real(obj.interface_heat_flux);
            obj.plant_absorbed_ir_sky = real(obj.plant_absorbed_ir_sky);
            obj.transpiration = real(obj.transpiration);
            obj.plant_convection = real(obj.plant_convection);
            obj.Qir = real(obj.Qir);
            obj.substrate_infrared_radiation = real(obj.substrate_infrared_radiation);
            obj.substrate_convection = real(obj.substrate_convection);
            obj.evaporation = real(obj.evaporation);
            obj.substrate_conduction = real(obj.substrate_conduction);
            obj.heating_load = real(obj.heating_load);
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