function [rho,Tinf,MM,N2,O2,OO,Ar,He,H,total_mass]=AirDenTemp(re,h,Date,Time,thick)

%% Air Density and Temperature 
% Calculate Latitude and Longitude
% Convert re (in ECI) to ECEF

% Make sure to check http://sail.msfc.nasa.gov/ for newest
% atmospheric model predictions. Download table 3, delete
% headers, and change file name in line 48

% Inputs: re = radius of orbit in ECI coordinates (m)
%         h = altitude (m)
%         Date = [day, month, year]; Use month number, name
%         Time = [hour, minute, second]; UTC, military time
%         thick = thickness flag, 1=high (95%), 2=medium (50%), 3=low (5%)
% Outpus: Tinf = K, Freestream static temp
%         rho = kg/m^3, air density
%         MM = Molecular mass
% Programs Called:
%   MET.m = translated from NASA fortran code for Jacchia atmosphere
%   ECI_to_ECEF.m = converts from ECI to ECEF; Included.

iday=[31,28,31,30,31,30,31,31,30,31,30,31];
year=365.2422;
%% Find Atmospheric Values for MET Input

iyr=Date(3);
mn=Date(2);
ida=Date(1);
ihr=Time(1);
min=Time(2);
sec=Time(3);

%Calculate decimal date
yr=iyr;
if mod(iyr,4) == 0
    if mod(iyr,100) ~= 0
        iday(2)=29; % Century not a leap year
    end
else
    iday(2)=28;
end

id=0;
if mn>1
    for iter=1:mn-1
        id=id+iday(iter);
    end
end
dd=id+ida;
dy=dd/year+yr;

% Look up f10 number and geomagnetic activity index
% Solar flux values measured in solar flux units (sfu): 1 sfu = 104 Jy = 10e-22 W/m2/Hz = 10e-19 erg/s/cm2/Hz
solar_flux=fopen('MSFC_Solar_Flux_Data_02_2019.txt');
data=fscanf(solar_flux,'%f %*s %f %f %f %f %f %f',[7 Inf]);
data=data';
fclose(solar_flux);

% Match data to correct row in data
Index=find(data(:,1)>dy,1,'first');
% Assign values based on percentiles
switch thick
    case 1 % High values, 95%
        f10=data(Index,2);
        f10b=mean(data(Index-2:Index+2,2));
        gi=data(Index,5);
    case 2 % Medium values, 50%
        f10=data(Index,3);
        f10b=mean(data(Index-2:Index+2,3));
        gi=data(Index,6);
    case 3 % low values, 5%
        f10=data(Index,4);
        f10b=mean(data(Index-2:Index+2,4));
        gi=data(Index,7);
end

%% Algorithm using MET
i1 = 2.0;
dcm = ECI_to_ECEF(iyr,mn,ida,ihr,min,sec);
recef=dcm*re;       % Earth radius in ECEF coords

lat=asin(recef(3)/norm(recef));
an_long = recef(1)/(norm(recef)*cos(lat));

if an_long > 1
    an_long = 1;
end

long = acos(an_long);
if isreal(long) == 0
    error
end

% Marshall Thermospheric Model (translated by Alexandra Long)
[~,Tinf,N2,O2,OO,Ar,He,H,MM,rho,~,~] = MET(h,lat,long,iyr,mn,ida,ihr,min,i1,f10,f10b,gi);
total_mass = N2+O2+OO+Ar+He+H;
N2 = N2/total_mass;
O2 = O2/total_mass;
OO = OO/total_mass;
Ar = Ar/total_mass;
He = He/total_mass;
H = H/total_mass;

end

function dcm = ECI_to_ECEF(yr,mon,da,hr,min,sec)
%% ECI to ECEF Direction Cosine Matrix Caculator
% This program calculates the transformation matrix from ECI to ECEF
% This uses the equations from "Aerospace Engineering: Orbital Mechanics
% for Engineering Students" by Howard Curtis
%
% Inputs: yr = year (1901-2099)
%         mon = month (1-12)
%         da = day (1-31)
%         hr = hour
%         min = minute
%         sec = seconds

%% Calculations
% Calculate Julian Date - Eq 5.48
J0 = 367*yr - fix(7*(yr + fix((mon + 9)/12))/4) ...
    + fix(275*mon/9) + da + 1721013.5;                  % J0 = Julian day number at 0 hr UT [days]

% Calculate T0 (Eq 5.49)
T0 = (J0 - 2451545)/36525;                              % T0 = time in Julian centuries between J0 and J2000 [-]

% Calculate thG0 (Eq 5.50)
thG0 = 100.4606184 + 36000.77004*T0 + 0.000387933*T0^2 - (2.583e-8)*T0^3; % Greenwich sidereal time at 0 hr UT [degrees]

% Make sure thG0 is between 0 and 360 degrees
if (thG0 >= 360)
    x = thG0 - fix(thG0/360)*360;
elseif (thG0 < 0)
    x = thG0 - (fix(thG0/360) - 1)*360;
end
thG0 = x; % degrees

% Calculate thG (Eq 5.51)
UT = hr + min/60 + sec/3600;            % Universal time [hrs]
thG = thG0 + 360.98564724*UT/24;        % Greenwich sidereal time at any UT [degrees]
thGr = deg2rad(thG);                    % convert to radians

% Direction cosine matrix - simple rotation about z axis
dcm = [ cos(thGr), sin(thGr),0;...
       -sin(thGr), cos(thGr),0;...
                0,         0,1];
end
