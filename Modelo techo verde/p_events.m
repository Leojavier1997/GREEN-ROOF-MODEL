
%EVENTS DURATION
%Efective precipitation events are separated by a certain number of hours

function [time_events]=p_events(Ptotal,t_s,dt)

%INPUT
%Ptotal         Total precipitacion (mm)
%t_s            Time interval between events (h)

%VARIABLE INITIATION 
%time_events indicates the amount of dt from the beginning of the events
time_events=zeros(length(Ptotal),1);

%The first value of the precipitation series starts the first event
for i=1:min(t_s/dt,length(Ptotal))
    time_events(i)=i;
end

%When the time between two records of total precipitation is greater than
%t_s, then that events ends and another begins
counter=t_s/dt;
for i=counter+1:length(Ptotal)
    if Ptotal(i)>0
        events0=0;
        k=1;
        for j=1:counter
            if Ptotal(i-k)>0
                events0=events0+1;
            end
            k=k+1;
        end
        if events0>0
            time_events(i)=time_events(i-1)+1;
        else
            time_events(i)=1;
        end
    else
        time_events(i)=time_events(i-1)+1;
    end
end


return





