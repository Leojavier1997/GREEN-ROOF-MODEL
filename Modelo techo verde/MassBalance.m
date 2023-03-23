%HYDROLOGICAL PROCESSES FOR PERVIOUS subarea

function [theta,DthetaDt,f,pe,red,F,Ft,Peffect,AWI,esc,irr_vol,Ptot_cum_event]= MassBalance(tstep,P,R,runon, E_T,general_inf, outflow_data, sub_data,plant_data2,theta,DthetaDt,f,pe,red,F,Ft,Peffect,AWI,esc,irr_vol,Ptot_cum_event)



%INPUT
%P                                      %Precipitation (mm)
%runon                                  %Runon volume coming from other subareas (mm)
%drainage                               %Drainage connection between this subarea and others (%)
%ET0                                    %Evapotranspiration (mm/h)
%R                                      %Irrigation (mm)
%time_inf                               %Temporal information
%input_sub                              %subarea parameters
%input_lay                              %All layers parameters
%Veg_data                             %Kcrop information for the vegetation involved


%Temporal information
dt=general_inf(2)/60;                         %Modeling time interval (h)
t_s=general_inf(3);                        %Time interval between events (h)

%Parameters for estimating hydrological processes

A=outflow_data(3);                         %subarea area (m2)
D=outflow_data(5);                         %Surface storage depth (mm)                                       
                     
nL=general_inf(4);                       %Number of layers in the subarea

vegcov=plant_data2(1)/100;          %Vegetation coverage
S=plant_data2(3);                   %Maximum water stored (mm)
k=plant_data2(2);                   %Leaf area index


%Layer's information
%Information for each layer is ordered in 'layer' matrix. Each columnn correspond to one layer
layer=nan(10,nL);
for j=1:nL
    layer(:,j)=sub_data(1:10);
end
%The parameters below are defined for all layers
d=layer(1,:)*1000/nL;                      %Layer thickness (mm)
theta_i=layer(2,:);                     %Initial soil water content  (m3/m3)
theta_r=layer(3,:);                     %Residual soil water content (m3/m3)
theta_s=layer(4,:);                     %Saturated soil water content (m3/m3)
n=layer(5,:);                           %Curve shape parameter
m=1-1./n;                               %Second curve shape parameter
L=layer(6,:);                           %Empirical pore tortuosity
Ks=layer(7,:);                          %Saturated hydraulic conductivity (mm/h)
psi_b=layer(8,:);                       %Bubbling pressure (mm)
%WP=layer(9,:);                          %Wilting point (m3/m3)
FC=layer(10,:);                         %Field capacity (m3/m3)
theta_e=theta_s(1)-theta_r(1);          %Effective soil water content (m3/m3)


%OUTPUT
%runoff                                 %Runoff volume that drains from this subarea (m3)
%flp                                    %Summary of the firts layer processes
    %Peffect                            %Precipitation plus runoff minus interception (mm)
    %F                                  %Cumulative infiltration (mm)
    %esc                                %effective surface runoff (mm)
    %f                                  %Infiltration rate (mm/h)
    %E_T                                %Evaporation or Evapotranspiration (mm/h)
%theta                                  %Soil water content (m3/m3)
%pe                                     %Percolation (mm/h)
%red                                    %Redistribution (mm/h)
%Q_sur                                  %Flow of surface runoff(m3/h)
%Q_sub                                  %Flow of subsurface runoff (m3/h)
%irr_vol                                %Irrigation volume (m3)

%VARIABLES'S INITIATION
e=zeros(length(P),1);                   %Variable assistant to calculate surface runoff (mm)
theta_p=zeros(1,nL);                    %Variable assistant to calculate soil water content (m3/m3)
K_p=zeros(nL,1);                        %Variable assistant to calculate unsaturated hydraulic conductivity (mm/h)
psi_p=zeros(nL,1);                      %Variable assistant to calculate suction head at wetting front (mm)

Aeq=zeros(nL,2+2*nL);                   %Matix assistant to calculate mass balance

for j=1:nL
    
    if j==1 %For the first layer
        Aeq(1,1)=dt/d(1);
        Aeq(1,2)=-dt/d(1);
        Aeq(j,3)=-dt/d(1);
        Aeq(1,3+nL)=-1;    
    else %For the deeper layers
        Aeq(j,j+2)=-dt/d(j);
        Aeq(j,j+1)=dt/d(j);
        Aeq(j,2+nL+j)=-1;
    end  
end

%TOTAL PRECIPITATION
Ptotal=P+runon/A*1000; %(mm)
%Ptotal includes rainfall plus runoff from upstream subarea
%Ptotal is uniformly distribuited


%INITIAL CONDITIONS
theta(1,:)=theta_i(:);                  %Initial soil water content (m3/m3)
i=tstep;

%IRRIGATION
Ri=R(i); 
irr_vol=irr_vol+Ri*A/1000; %Irrigation volume (m3)
    
%INTERCEPTION (mm) 
    %Total precipitation minus interception
        Peffect(i)=Ptotal(i)+Ri;
        
    %What are the new precipitation events?
        %time_events indicates the amount of dt from the beginning of the events
        [time_events]=p_events(Peffect(1:i),t_s,dt);
        %Interception
        if i==1 %At the beginning
            ET0_p=E_T(i);
            time_events_p=time_events(i);
            Ptot_cum_event=Peffect(i);
        else
            ET0_p=[E_T(i) E_T(i-1)];
            time_events_p=[time_events(i) time_events(i-1)];    
        end 
        [Peffect(i),Ptot_cum_event]=Interception(Peffect(i),ET0_p,Ptot_cum_event,time_events_p,k,S*vegcov,dt);
           
       
%AVAILABLE WATER TO INFILTRATE (mm)
%If there is a reservoir, AWI includes the amount of water stored in the previous period
    AWI(i)=AWI(i)+Peffect(i);
    
   
%PERCOLATION (mm/h) - For all layers    
    for j=1:nL
        %van Genuchten parameters
        Se=(theta(i,j)-theta_r(j))/(theta_s(j)-theta_r(j));     %Relative saturation(mm3/mm3)
        K=(Ks(j)*Se^L(j))*(1-(1-Se^(1/m(j)))^(m(j)))^2;         %Unsaturated hydraulic conductivity (mm/h)
        t=max((theta(i,j)-FC(j))*d(j)/K,0);                     %Travel time through the layer(h)                        
        %Percolation (mm/h)
        if theta(i,j)>FC(j)
            pe(i,j)=real((theta(i,j)-FC(j))*(1-exp(-dt/t))*d(j)/dt); %se multipica por d/dt para tener consistencia en las unidades
        else
            pe(i,j)=0;
        end
    end

%INFILTRATION (mm/h) - for the first layer
    %Updating new precipitation events
    [time_events]=p_events(AWI(1:i),t_s,dt);
    F00=0;
    %Initial flooding for Green Ampt calculation
    if i==1
         F0=F00;
    %No flooding is assumed for each new event
    elseif time_events(i)==1 
         F0=0;
    else
    %Initial flooding is equal to cumulative infiltration from previous time interval
         F0=Ft(i-1);
    end
    
    %Suction head at wetting front (mm)
    psi=psi_b(1)*(((theta(i,1)-theta_r(1))/(theta_s(1)-theta_r(1)))^(-1/m(1))-1)^(1/n(1));
    %Infiltration (mm/h)  
   
    [f(i),Ft(i)]=Green_Ampt(F0,psi,Ks(1),dt,Se(1),theta_e,AWI(i));

    
%MASS BALANCE - For all layers
    limit_sup=0; %Auxiliary variable
    for j=1:nL
        %Soil water content considering the above hydrological processes
        if j==1 %First layer
            theta_p(j)=theta(i,j)+(f(i)-E_T(i)-pe(i,j))*dt/d(j);
        else %Deeper layers
            theta_p(j)=theta(i,j)+(pe(i,j-1)-pe(i,j))*dt/d(j);
        end
        %If the soil water content of the first layer exceeds theta_s, then
        %the infiltration rate is decreased
        if theta_p(1)>theta_s(1)
            Ft(i)=Ft(i)-(f(i)-((theta_s(1)-theta(i,1))*d(1)/dt+E_T(i)+pe(i,1)))*dt;
            f(i)=((theta_s(1)-theta(i,1))*d(1)/dt+E_T(i)+pe(i,1));
        end
        %Try again
        if j==1 %First layer
            theta_p(j)=theta(i,j)+(f(i)-E_T(i)-pe(i,j))*dt/d(j);
        else %Deeper layers
            theta_p(j)=theta(i,j)+(pe(i,j-1)-pe(i,j))*dt/d(j);
        end
        %Analysis of the limits of soil water content: 
        %Register if theta of some of the layers exceeds theta_s or is less than theta_r
        if theta_p(j)<theta_r(j) || theta_p(j)>theta_s(j)
            limit_sup=1;
        end
    end
    %If theta is maintained within the limits, its value is updated
    if limit_sup==0 
        theta(i+1,:)=theta_p;
    %If theta exceeds a limit, hydrological rates are adjusted as optimization problem     
    else
        %Initial values
        x0=[f(i),E_T(i),pe(i,:),theta_p];   %x=[f E_T pe theta]
        %Limits
        li=[0,0,zeros(1,nL),theta_r*1.01]; %Lower limit
        ls=[f(i),E_T(i),pe(i,:),theta_s*0.99]; %Upper limit
        options=optimset('Algorithm','active-set');
        %The problem is solved by maximizing the sum of rates f+E_T+pe
        [x]=fmincon(@(x)(-sum(x(1:2+nL))),x0,[],[],Aeq,-theta(i,:),li,ls,[],options);        
        %Rates are updated
        %If theta exceeds theta_s, the infiltration rate must be reduced
        if x(1)<f(i)
            Ft(i)=Ft(i)-(f(i)-x(1))*dt;
            f(i)=x(1);
        end
        %If theta is less than theta_r, the evapotranspiration or
        %percolation rate must be reduced
        %E_T(i)=x(2);
        pe(i,:)=x(3:2+nL);
        %The rate of change of soil water content is updated
        for j=1:nL
            if j==1 %First layer
                theta(i+1,j)=theta(i,j)+(f(i)-E_T(i)-pe(i,j))*dt/d(j);
            else %Deeper layers
                theta(i+1,j)=theta(i,j)+(pe(i,j-1)-pe(i,j))*dt/d(j);
            end
        end
    end
    
    
%CUMULATIVE INFILTRATION (mm) - For the first layer
%Infiltration depth during the time step dt. It value is calculated based on
%cumulative infiltration for each precipitation events
    if i==1
        dF=Ft(i)-F00;
    elseif time_events(i)==1 %At the beginning of the event
        dF=Ft(i);
    else %During the event
        dF=Ft(i)-Ft(i-1);
    end
    %Cumulative infiltration considering the complete simulation
    if i~=1
        F(i)=real(F(i-1)+dF);
    else
        F(i)=Ft(i);
    end
    
    
%SURFACE RUNOFF (mm) - For the first layer
    %Infiltration depth during the time step dt. It value is calculated based on
    %cumulative infiltration considering the complete simulation
    if i==1 %At the beginning of the simulation
        Fp=F(i);
    else %During the simulation
        Fp=F(i)-F(i-1);
    end
    %The rainfall excess (mm) is calculated in each time as the difference
    %between AWI and the infiltration depth
    e(i)=max(AWI(i)-Fp,0);
    %If the subarea no considers a surface storage, the runoff is
    %equal to rainfall excess
    if D==0
        esc(i)=e(i);
    %If the subarea considers a surface storage, then surface runoff
    %is generated only once the storage capacity is full
    elseif D>0 
       if e(i)>D %If storage is full
           esc(i)=e(i)-D;
           AWI(i+1)=D;
       else %If storage is not full, no runoff is generated
           AWI(i+1)=e(i);
           esc(i)=0;
       end
    end
    
%REDISTRIBUTION (mm/h) - for all layers
    %Redistribution occurs when the infiltration is null and the
    %subarea is composed of more than one layer
    if f(i)<1e-12 && nL>1
        rr=ones(1,nL-1); %Redistribution factor
        for j=1:nL
            %van Genuchten parameters considering the theta values after
            %considering the above hydrological processes
            Se=(theta(i+1,j)-theta_r(j))/(theta_s(j)-theta_r(j));
            psi_p(j)=real(psi_b(j)*(Se^(-m(j)^(-1))-1)^(n(j)^(-1)));
            K_p(j)=real(Ks(j)*Se^L(j)*(1-(1-Se^(n(j)/(n(j)-1)))^((n(j)-1)/n(j)))^2);
        end
        for j=2:nL-1    
            %If the soil is composed of 3 layers or more and the flow from the middle layers is
            %established both upward and downward, the total redistribution
            %flow obtained is split according to a factor 'rr'
            if psi_p(j-1)>psi_p(j) && psi_p(j+1)>psi_p(j)
                rr(j-1)=abs(psi_p(j)-psi_p(j-1))/(abs(psi_p(j)-psi_p(j-1))+abs(psi_p(j)-psi_p(j+1)));
                rr(j)=abs(psi_p(j)-psi_p(j+1))/(abs(psi_p(j)-psi_p(j-1))+abs(psi_p(j)-psi_p(j+1)));
            end
        end
        for j=1:nL-1
            if psi_p(j)<=psi_p(j+1) %The flow downward
                red(i,j)=K_p(j)*((psi_p(j)-psi_p(j+1))/(0.5*(d(j)+d(j+1)))-1)*rr(j);%(negative)
            else %The flow upward
                red(i,j)=K_p(j+1)*((psi_p(j)-psi_p(j+1))/(0.5*(d(j)+d(j+1)))-1)*rr(j);%(positive)
            end
        end
        limit_sup=0; %Auxiliary variable
        for j=1:nL
            if j==1 %For the firt layer
                theta_p(j)=theta(i+1,j)+red(i,j)*dt/d(j);
            elseif j>1 && j<nL %For the middle layers
                theta_p(j)=theta(i+1,j)+(red(i,j)-red(i,j-1))*dt/d(j);
            elseif j==nL %For the last layer
                theta_p(j)=theta(i+1,j)-red(i,j-1)*dt/d(j);
            end
            %Analysis of the limits of soil water content: 
            %Register if theta of some of the layers exceeds theta_s or is less than theta_r
            if theta_p(j)<theta_r(j) || theta_p(j)>theta_s(j)
                limit_sup=1;
            end
        end
        %If theta exceeds a limit, redistribution rate is adjusted as
        %optimization problem  
        if limit_sup==1  
            dir_red=sign(red(i,1:nL-1)); %Redistribution flow orientation
            %Matrix Aer
            Aer=zeros(nL,2*nL-1);
            for j=1:nL
                if j==1 %For the first layer
                    Aer(j,j)=1;
                    Aer(j,nL+j)=-dir_red(j)*dt/d(j); 
                elseif j>1 && j<nL %For the middle layers
                    Aer(j,j)=1;
                    Aer(j,nL+j-1)=dir_red(j-1)*dt/d(j);
                    Aer(j,nL+j)=-dir_red(j)*dt/d(j);
                else %For the last layer
                    Aer(j,j)=1;
                    Aer(j,nL+j-1)=dir_red(j-1)*dt/d(j);
                end
            end
            %Initial values
            y0=[theta_p,abs(red(i,1:nL-1))];   %y=[theta(nL) red(nL-1)]
            %Limits
            li=[theta_r*1.01,zeros(1,nL-1)]; %Lower limit
            ls=[theta_s*0.99,abs(red(i,1:nL-1))]; %Upper limit
            options=optimset('Algorithm','active-set');
            %The problem is solved by maximizing the sum of rates red
            [y]=fmincon(@(y)(-sum(y(nL+1:end).*(1-dir_red)/2)),y0,[],[],Aer,theta(i+1,:),li,ls,[],options);
            %Rates are updated
            red(i,1:nL-1)=dir_red.*y(nL+1:end);               
        end    
    end
  
%RATE OF CHANGE OF SOIL WATER CONTENT
    for j=1:nL
        if j==1
            DthetaDt(i,j)=f(i)-E_T(i)-pe(i,j)+red(i,1);
        else
            DthetaDt(i,j)=pe(i,j-1)-pe(i,j)-red(i,j-1)+red(i,j);
        end
        theta(i+1,j)=theta(i,j)+DthetaDt(i,j)*dt/d(j);
    end
    

return

