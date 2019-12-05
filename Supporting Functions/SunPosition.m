function [R_sun] = SunPosition(Epoch)
%% SunPosition.m

%Purpose: Determine the location of the Sun relative to the Earth in ECI
%using an analytical approximation from the Astronomical Almanac (1992)

%Created:  Sebastian Tamrazian 1/7/2018
 
%Inputs:    Epoch - %mm-dd-yyyy HH:MM:SS' [string format]
%Outputs:   R_sun - sun vector from earth in ECI

%Variables: 
%   Year,Month,Day,Hour,Min,Sec - doubles that describe the date
%   Juliandate - Julian date at the current epoch

%Supporting Functions: 
%   None

[Year,Month,Day,Hour,Min,Sec] = datevec(datenum(Epoch));

Juliandate = 367*Year - floor((7*(Year + floor((Month + 9)/12)))/4) + ...
    floor((275*Month)/9) + Day + 1721013.5 + (Hour + Min/60 + Sec/3600)/24;

Juliancenturies = (Juliandate-2451545)/36525;

Mean_longitude = 280.46+36000.771*Juliancenturies; %deg

Mean_anomaly = 357.5291092 +35999.05034*Juliancenturies; %deg

Ecliptic_longitude = Mean_longitude + 1.914666471*sind(Mean_anomaly)+.019994643*sind(2*Mean_anomaly); %deg

Distance_sun = 1.000140612-.016708617*cosd(Mean_anomaly)-.000139589*cosd(2*Mean_anomaly); %AU

Ecliptic_obliquity = 23.439291-.0130042*Juliancenturies; %deg

Distance_sun = Distance_sun * 149597871e3; %convert from AU to m

R_sun = [Distance_sun*cosd(Ecliptic_longitude);...
         Distance_sun*cosd(Ecliptic_obliquity)*sind(Ecliptic_longitude);...
         Distance_sun*sind(Ecliptic_obliquity)*sind(Ecliptic_longitude)];
end

