% 

clear all 
close all
clc

%%%             LLENAR ESTO             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %input_file="INPUTSLive_AgostoC.xlsx";
    %input_file="INPUTSLive_AgostoA.xlsx";
    %input_file="INPUTSMelbourne10cm_revgf.xlsx";
    input_file="INPUTSChicago_AR_Prueba2.xlsx";

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



surface_sensor_depth = 0.000;
mid_sensor = 0.05;
deep_sensor = 0.1;


%% PARSE GENERAL DATA

% Read
[numbers,files] = xlsread(input_file,'General','K4:K14');

% Parse Files
out_file = files(1);
if(files(2) ~= "")    
    use_epw = true;    
    epw_file = files(2);
    disp("EPW file '" +epw_file+"' will be used");
else
    use_epw = false;
end
model_to_use = files(4);
location = files(7);
clear files

% parse numbers
meteo_dt = numbers(1); %minutes
sim_dt = numbers(2); %minutes
event_dt = numbers(3); %hours
n_layers = numbers(4);
max_steps = numbers(7);
Area = numbers(9);
roof_height = numbers(10);

% These values correspond to OCEAN

if location == "Country"
    SiteWindBLHeight = 270;
    SiteWindExp = 0.14;    
elseif location == "Suburbs"
    SiteWindBLHeight = 370;
    SiteWindExp = 0.22;
elseif location == "City"
    SiteWindBLHeight = 460;
    SiteWindExp = 0.33;    
elseif location == "Ocean"    
    SiteWindBLHeight = 210;
    SiteWindExp = 0.1;    
elseif location == "Urban"
    SiteWindBLHeight = 370;
    SiteWindExp = 0.22;    
else
    error('Not recognized location');
end
clear numbers

%% PARSE OUTFLOW


%% PARSE OUTDOOR 
interior_temperature = xlsread(input_file,"TODO","A:A")+273.15;
measured_VWC = xlsread(input_file,"TODO","B:B");
if(use_epw)
    % Load data from EPW
    weatherdata = csvread(char(epw_file),8,6);
    weatherdata = [weatherdata(end,:);weatherdata];
else
    % Load data from excel
    weatherdata = xlsread(input_file,"Outdoors","F:K");
end


%% CREATE PLANT
plant = Plant();

plant_data = xlsread(input_file,'Plants','E5:E14');

plant.LAI=plant_data(1);        %leaf area index A FEB Tabares  
plant.rho=plant_data(3);     %sw reflectivity plant
plant.em=plant_data(4);    %emissivity planta
plant.k=plant_data(5);      % thermal conductivity planta
plant.height=plant_data(6);   %  height of plant [m]
plant.ks=plant_data(7); % Extinsion coefficient   
plant.ks_ir = plant_data(8); % IR extinsion coefficient
plant.rsmin=plant_data(9);      %[s/m] resistencia estomatica 
plant.z=2;       % [m] height of wind measurements
plant.Zog = plant_data(10); % VARIA SEGUN SMOOOTH
clear plant_data

%% CREATE SUBSTRATE
substrate = Substrate();
sub_data = xlsread(input_file,'Substrate','E4:E18');

substrate.depth=sub_data(1);     % depth substrate  (ex t)
initial_VWC = sub_data(2); % it is used later
substrate.VWCresidual = sub_data(3);
substrate.VWCsat=sub_data(4);     % VWC saturated
substrate.VWCwilting=sub_data(9);    % VWC wilting point
substrate.VWCfc = sub_data(10);
substrate.rho=sub_data(11);     %sw reflectivity substrate
substrate.k=sub_data(12); %W/mK;
substrate.dens =sub_data(13); %kg/m^3
substrate.Cp=sub_data(14); % J/kgK;
substrate.em = sub_data(15);
substrate.phi=substrate.VWCsat; 
% substrate.Lsubm=0.05;
clear sub_data


%% CREATE ROOF

roof = Roof();
roof_data = xlsread(input_file,'Support','E4:E7');

roof.depth = roof_data(1); % meters
roof.k = roof_data(2); %thermal conductivity
roof.density = roof_data(3); 
roof.Cp = roof_data(4); % heat capacity

clear roof_data

%% CREATE MODEL
disp("Using model "+model_to_use);
if model_to_use == "Tabares"
    model = TabaresThermalMass(plant,substrate,roof,n_layers,n_layers,Area,sim_dt);
elseif model_to_use == "Sailor"
    model = SailorThermalMass(plant,substrate,roof,n_layers,n_layers,Area,sim_dt);
else    
    error('Unkown model to use');
end
model.VWC = initial_VWC*ones(n_layers,1);

%% FIX TIMESTEPS
if use_epw 
    meteo_dt = 60; %EPW has data every 60 minutes
end


%% Calculate number of iterations
n_sub_tsteps = meteo_dt/sim_dt;
total_steps = n_sub_tsteps*(max_steps-1);


%% CREATE RESULTS VECTORS
result_T_sky = zeros(total_steps,1);
result_T_out = zeros(total_steps,1);
result_wind_speed = zeros(total_steps,1);
result_Rsh = zeros(total_steps,1);
result_T_plants = zeros(total_steps,1);
result_T_substrate = zeros(total_steps,1);
result_T_interior = zeros(total_steps,1);
result_VWC = zeros(total_steps,1);
result_interface_heat_flux = zeros(total_steps,1);
result_evaporation = zeros(total_steps,1);
result_transpiration = zeros(total_steps,1);
result_substrate_convection = zeros(total_steps,1);
result_plant_convection = zeros(total_steps,1);
result_Rain = zeros(total_steps,1);
result_ET = zeros(total_steps,1);

result_rs = zeros(total_steps,1);

result_fsolar = zeros(total_steps,1);
result_fvpd = zeros(total_steps,1);
result_fvwc = zeros(total_steps,1);
result_ftemp = zeros(total_steps,1);

result_ra = zeros(total_steps,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_5cm = zeros(total_steps,1);
T_10cm = zeros(total_steps,1);
T_15cm = zeros(total_steps,1);
T_20cm = zeros(total_steps,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plants balance
result_plant_absorbed_solar = zeros(total_steps,1);
result_plant_absorbed_ir_sky = zeros(total_steps,1);      
transpiration = zeros(total_steps,1);

% Both
result_Qir = zeros(total_steps,1);

% Substrate balance
result_substrate_solar_radiation = zeros(total_steps,1);
result_substrate_infrared_radiation = zeros(total_steps,1);
result_substrate_conduction = zeros(total_steps,1);

%% SIMULATE
if use_epw    
    sim_desc = 'Performing simulation with inputs from '+input_file+' using '+model_to_use+' model and data from '+epw_file;
else    
    sim_desc = 'Performing simulation with inputs from '+input_file+' using '+model_to_use+' model and custom data';
end
h = waitbar(0,char(sim_desc));
for main_step = 1:length(weatherdata)
    
    this_data_line = weatherdata(main_step,:);
    next_data_line = weatherdata(main_step+1,:); 
    
    for k=1:n_sub_tsteps
                    

        tstep = (main_step-1)*n_sub_tsteps+k;
        waitbar( tstep/ (max_steps*n_sub_tsteps));
            
        % interpolate   
        this_inner_t = interior_temperature(main_step);
        next_inner_t = interior_temperature(main_step+1);
        inner_T = this_inner_t + (k-1)*(next_inner_t - this_inner_t)/n_sub_tsteps;        
        model.T_interior = inner_T;                                                
        
        if(use_epw)
            data_line = EPWWeatherDataLine(this_data_line + (k-1)*(next_data_line - this_data_line)/n_sub_tsteps,roof_height,SiteWindBLHeight,SiteWindExp);
        else
            data_line = WeatherDataLine(this_data_line + (k-1)*(next_data_line - this_data_line)/n_sub_tsteps);
        end
        
        % Advance
        model = model.moveForward(data_line);
        
        % Update moisture
        this_vwc = measured_VWC(main_step);
        next_vwc = measured_VWC(main_step+1);
        vwc_now = this_vwc+ (k-1)*(next_vwc- this_vwc)/n_sub_tsteps;                        
        model.VWC = vwc_now*ones(n_layers,1);        
        
        result_ET(tstep) = model.et_mm_hour;
        
        % Retrieve and store data
        result_T_sky(tstep) = data_line.Tsky;
        result_T_out(tstep) = data_line.Tair;
        result_wind_speed(tstep) = data_line.U;
        result_Rsh(tstep) = data_line.R_sh;
        result_Rain(tstep) = data_line.rainfall;
        result_T_plants(tstep) = model.T_plants;
        result_T_substrate(tstep) = model.getTemperature(surface_sensor_depth);
        result_T_interior(tstep) = model.T_interior;
        result_VWC(tstep) = model.VWC(1);
        result_interface_heat_flux(tstep) = model.interface_heat_flux;
        result_evaporation(tstep) = model.evaporation;
        result_transpiration(tstep) = model.transpiration;
        result_substrate_convection(tstep) = model.substrate_convection;
        result_plant_convection(tstep) = model.plant_convection;
        
        result_rs(tstep) = model.rs;
        
        result_fsolar(tstep)=model.fsolar;
        result_fvpd(tstep) = model.fvpd;
        result_fvwc(tstep) = model.fvwc;
        result_ftemp(tstep) = model.ftemp;

        result_ra(tstep) = model.ra;

        
        result_plant_absorbed_solar(tstep) = model.plant_absorbed_solar ;
        result_plant_absorbed_ir_sky(tstep) = model.plant_absorbed_ir_sky;      
        result_Qir(tstep) = model.Qir;
        result_substrate_solar_radiation(tstep) = model.substrate_solar_radiation;
        result_substrate_infrared_radiation(tstep) = model.substrate_infrared_radiation;
        result_substrate_conduction(tstep) = model.substrate_conduction;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T_5cm(tstep) = model.getTemperature(mid_sensor);
        T_10cm(tstep) = model.getTemperature(deep_sensor);
        T_15cm(tstep) = model.getTemperature(0.15);
        T_20cm(tstep) = model.getTemperature(0.195);
        
        result_T_substrate(tstep) = model.getTemperature(0.005);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    % break if needed
    if main_step >= (max_steps-1)
        break
    end
end
close(h);

headers = [
    "Sky temperature (°C)",...
    "Exterior temperature (°C)",...
    "Wind speed (m/s)",...
    "Global solar radiation (W/m2)",...
    "Substrate Surface Temperature (°C)",...
    "Foliage Temperature (°C)",...
    "VWC",...
    "VWC",...
    "Substrate sensible heat transfer (W/m2)",...
    "Foliage sensible heat transfer (W/m2)",...
    "Foliage Transpiration (W/m2)",...
    "Substrate evaporation (W/m2)",...
    "Rainfall (mm)",...
    "Interior temperature (°C)",...
    "Heat flux below surface (W/m2)",...
    "Evapotranspiration (mm/hour)",...
    "plant_absorbed_solar",...
    "result_plant_absorbed_ir_sky",...     
    "result_Qir",...
    "result_substrate_solar_radiation",...
    "result_substrate_infrared_radiation",...
    "result_substrate_conduction",...
    "ra",...
    "fsolar",...
    "fvpd",...
    "fvwc",...
    "ftemp",...
    "rs"...
    ];

results = [...
    result_T_sky-273.15,...
    result_T_out-273.15,...
    result_wind_speed,...
    result_Rsh,...
    result_T_substrate-273.15,...
    result_T_plants-273.15,...
    result_VWC,...
    result_VWC,...
    result_substrate_convection, ...
    result_plant_convection, ...
    result_transpiration, ...
    result_evaporation, ...
    result_Rain, ...
    result_T_interior-273.15,...
    result_interface_heat_flux,...
    result_ET,...    
    result_plant_absorbed_solar,...
    result_plant_absorbed_ir_sky,...     
    result_Qir,...
    result_substrate_solar_radiation,...
    result_substrate_infrared_radiation,...
    result_substrate_conduction,...
    result_ra,...
    result_fsolar,...
    result_fvpd,...
    result_fvwc,...
    result_ftemp,...
    result_rs...

];


xlswrite(char(out_file)+".xlsx",[headers;real(results(1:n_sub_tsteps:end,:))],char(model_to_use));

headers2 = ["5cm","10cm","15cm","20cm"];
deep_results = [T_5cm(1:n_sub_tsteps:end,:), T_10cm(1:n_sub_tsteps:end,:), T_15cm(1:n_sub_tsteps:end,:), T_20cm(1:n_sub_tsteps:end,:)];
xlswrite(char(out_file)+".xlsx",[headers2; real(deep_results)-273.15],char(model_to_use)+"_deep");






