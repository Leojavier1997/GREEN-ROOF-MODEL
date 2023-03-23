% 

clear all 
close all
clc

%%%             LLENAR ESTO             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %input_file = "INPUTSChicago_short.xlsx";
    %input_file = "INPUTSLIVE_v03.xlsx";
    %input_file = "INPUTSLIVE_v07.xlsx";
    input_file = "INPUTSLive_AgostoC.xlsx";
    surface_sensor_depth = 0.03;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% GF
general_inf=numbers(1:4);
nL=n_layers; 
clear numbers

%GF
Dt=meteo_dt/60; % Dt expressed in hours                      
dt=sim_dt/60; % Dt expressed in hours  

%% PARSE OUTFLOW
%GF
outflow_data=xlsread(input_file,'Outflow','K4:K11');

%% PARSE OUTDOOR 
interior_temperature = xlsread(input_file,"TODO","A:A")+273.15;

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

%GF
plant_data2=xlsread(input_file,'Plants','E4:E14');

%% CREATE SUBSTRATE
substrate = Substrate();
sub_data = xlsread(input_file,'Substrate','E4:E18');

substrate.depth=sub_data(1);     % depth substrate  (ex t)
initial_VWC = sub_data(2); % This is used later.
substrate.VWCresidual = sub_data(3);
substrate.VWCsat=sub_data(4);     % VWC saturated
substrate.VWCwilting=sub_data(9);    % VWC wilting point
substrate.VWCfc = sub_data(10);
substrate.rho=sub_data(11);     %sw reflectivity substrate
substrate.k=sub_data(12); %W/mK;
substrate.dens =sub_data(13); %kg/m^3
substrate.Cp=sub_data(14); % J/kgK;
substrate.em = sub_data(15);
% substrate.phi=0.85; 
% substrate.Lsubm=0.05;
%clear sub_data
%GF : GM--> Initial VWC is given in cell E5
%GM: De acuerdo... modificado 
%substrate.VWC = 0.23; % initial Volumetric Water Content


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
n_sub_tsteps = meteo_dt/sim_dt;
total_steps = n_sub_tsteps*(max_steps-1);

%% INTERPOLATION FOR PRECIPITATION P AND IRRIGATION R
%GF
for i=1:length(weatherdata(:,1))
    if(use_epw)
        P(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherdata(i,28)*dt/Dt;
        R(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherdata(i,6)*0;
    else
        P(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherdata(i,1)*dt/Dt;
        R(1+round((i-1)*Dt/dt):round(i*Dt/dt),1)=weatherdata(i,6)*dt/Dt;
    end
end

%% CREATE RESULTS VECTORS
result_T_sky = zeros(total_steps,1);
result_T_out = zeros(total_steps,1);
result_wind_speed = zeros(total_steps,1);
result_Rsh = zeros(total_steps,1);
result_T_plants = zeros(total_steps,1);
result_T_substrate = zeros(total_steps,1);
result_T_5cm = zeros(total_steps,1); %5cm
result_T_10cm = zeros(total_steps,1); %10cm
result_T_15cm = zeros(total_steps,1); %15cm
result_T_interior = zeros(total_steps,1);
result_VWC_surface = zeros(total_steps,1);
result_VWC_mid = zeros(total_steps,1);
result_interface_heat_flux = zeros(total_steps,1);
result_evaporation = zeros(total_steps,1);
result_transpiration = zeros(total_steps,1);
result_substrate_convection = zeros(total_steps,1);
result_plant_convection = zeros(total_steps,1);
result_Rain = zeros(total_steps,1);
result_ET = zeros(total_steps,1);
result_heating_load = zeros(total_steps,1);

% Plants balance
result_plant_absorbed_solar = zeros(total_steps,1);
result_plant_absorbed_ir_sky = zeros(total_steps,1);      

% Both
result_Qir = zeros(total_steps,1);

% Substrate balance
result_substrate_solar_radiation = zeros(total_steps,1);
result_substrate_infrared_radiation = zeros(total_steps,1);
result_substrate_conduction = zeros(total_steps,1);

%% INITIALISATION OF VARIABLES
%GF
runon=zeros(length(P),1);               %Runoff volume entering to the subcatchment (m3)
runoff=zeros(length(P),1);              %Runoff volume which drains from the subcatchment (m3)
Q_out=zeros(fix(length(P)*1.1),1);      %Outflow from the subcatchment to the street (m3/h) 

theta=zeros(length(P),nL);              %Soil water content (m3/m3)
DthetaDt=zeros(length(P),nL);           %Rate of change in soil water content (mm/h)
f=zeros(length(P),1);                   %Infiltration (mm/h)
pe=zeros(length(P),nL);                 %Percolation (mm/h)

red=zeros(length(P),nL);                %Redistribution (mm/h)
F=zeros(length(P),1);                   %Cumulative infiltration (mm)
Ft=zeros(length(P),1);                  %Cumulative infiltration to calculate Green Ampt (mm)
Peffect=zeros(length(P),1);             %Precipitation plus runoff minus interception (mm)
AWI=zeros(length(P),1);                 %Available water to infiltrate (mm)
esc=zeros(length(P),1);                 %effective surface runoff (mm)

%% INITIAL VALUES
%GF
irr_vol=0;                              %Irrigated volume (m3)
Ptot_cum_event=0;

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
                        
        % GF : to be completed by GM
        result_ET(tstep) = model.et_mm_hour; % Evaporation + Transpiration in mm/hour
        %E_T=... %GF expressed in mm/h. Must be a matrix (1 line per timestep)
    
        % Perform the mass balance and determine the VWC (Volumetric Water
        % Content)of each layer
 
        [theta,DthetaDt,f,pe,red,F,Ft,Peffect,AWI,esc,irr_vol,Ptot_cum_event]...
            = MassBalance(tstep,P,R,runon, result_ET,general_inf, outflow_data,... 
                sub_data,plant_data2,theta,DthetaDt,f,pe,red,F,Ft,Peffect, AWI,...
                esc,irr_vol,Ptot_cum_event);
       
        % GF : to be modified by GM 
        model.VWC= theta(tstep,:);            
        
        
        % Retrieve and store data
        result_T_sky(tstep) = data_line.Tsky;
        result_T_out(tstep) = data_line.Tair;
        result_wind_speed(tstep) = data_line.U;
        result_Rsh(tstep) = data_line.R_sh;
        result_Rain(tstep) = data_line.rainfall;
        result_T_plants(tstep) = model.T_plants;
        result_T_substrate(tstep) = model.getTemperature(surface_sensor_depth);
        result_T_5cm(tstep) = model.getTemperature(0.05); %5cm
        result_T_10cm(tstep) = model.getTemperature(0.1); %10cm
        result_T_15cm(tstep) = model.getTemperature(0.15); %15cm
        result_T_interior(tstep) = model.T_interior;
        result_VWC_surface(tstep) = model.VWC(1);
        result_VWC_mid(tstep) = model.getVWC(substrate.depth/2);
        result_interface_heat_flux(tstep) = model.interface_heat_flux;
        result_evaporation(tstep) = model.evaporation;
        result_transpiration(tstep) = model.transpiration;
        result_substrate_convection(tstep) = model.substrate_convection;
        result_plant_convection(tstep) = model.plant_convection;
        result_heating_load(tstep) = model.heating_load;
        
        result_plant_absorbed_solar(tstep) = model.plant_absorbed_solar ;
        result_plant_absorbed_ir_sky(tstep) = model.plant_absorbed_ir_sky;      
        result_Qir(tstep) = model.Qir;
        result_substrate_solar_radiation(tstep) = model.substrate_solar_radiation;
        result_substrate_infrared_radiation(tstep) = model.substrate_infrared_radiation;
        result_substrate_conduction(tstep) = model.substrate_conduction;
        
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
    "Temperature 5cm deep (°C)",...
    "Temperature 10cm deep (°C)",...
    "Temperature 15cm deep (°C)",...
    "Foliage Temperature (°C)",...
    "VWC surface",...
    "VWC mid depth",...
    "Substrate sensible heat transfer (W/m2)",...
    "Foliage sensible heat transfer (W/m2)",...
    "Foliage Transpiration (W/m2)",...
    "Substrate evaporation (W/m2)",...
    "Rainfall (mm)",...
    "Interior temperature (°C)",...
    "Heat flux below substrate (W/m2)",...
    "Evapotranspiration (mm/hour)",...
    "plant_absorbed_solar",...
    "plant_absorbed_ir_sky",...     
    "Qir",...
    "substrate_solar_radiation",...
    "substrate_infrared_radiation",...
    "substrate_conduction",...
    "heating_load (W/m2)"
    ];

results = [...
    result_T_sky-273.15,...
    result_T_out-273.15,...
    result_wind_speed,...
    result_Rsh,...
    result_T_substrate-273.15,...
    result_T_5cm-273.15,...
    result_T_10cm-273.15,...
    result_T_15cm-273.15,...
    result_T_plants-273.15,...
    result_VWC_surface,...
    result_VWC_mid,...
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
    result_heating_load
];

xlswrite(char(out_file),[headers;results(1:n_sub_tsteps:end,:)],char(model_to_use));




