 function [ Freestream ] = AtmosphericProperties( EARTH, Time_epoch, Position, Velocity )
%% AtmosphericProperties.m function

%Purpose: Determine freestream properties surrounding spacecraft using
%Marshall Engineering Thermosphere model and Dr. Alexandra Long's
%translated code (in supporting functions)

%Created:  Sebastian Tamrazian 1/20/2019

% Inputs
%   EARTH: structure with fields:
%                     1. .SGP : Earth standard gravitational parameter [m^3/s^2]
%                     2. .MASS : Earth mass [kg]
%                     3. .EQRADIUS : Earth equatorial radius [m]
%                     4. .PORADIUS : Earth polar radius [m]
%                     5. .J2 : first zonal harmonic coefficient of Earth [-]
%
%   Time_epoch : time string in format 'dd-mm-yyyy HH:MM:SS'
%   Position : r vector [m]
%   Velocity : v vector [m/s]

% Outputs
%   Freestream
%             .density          -[kg/m^3]
%             .temperature      -[K]
%             .R                -[J/(mol K)]
%             .velocity         -[m/s]

% Supporting Functions: 
%    AirDenTemp.m

% Ideal gas constant
R_UNIVERSAL           = 8.3144598;       %[J/(mol K)]
rotrate               = 7.2921159E-5;    %[rad/s]
% rotrate = 0;

% function requires data and time %as numerical vector inputs to find solar 
% flux properties from indexed database from MSFC

[Datevector,Timevector] = strtok(Time_epoch); %Parse the cell

Datevector = str2double(strsplit(Datevector,'-')); % Parse the string
Timevector = str2double(strsplit(Timevector,':')); % Parse the string

thick = 2; % Values 1-3, 1 being thickest atmosphere

Altitude = norm(Position)-EARTH.EQRADIUS; %required distances in m

% Obtain Freestream conditions using Jacchia 70 model and Marshall
% Engineering Thermosphere model

[Freestream.density,Freestream.temperature,MM,~,~,~,~,~,~,~]=AirDenTemp(EARTH.EQRADIUS,Altitude,Datevector,Timevector,thick);

Freestream.R = R_UNIVERSAL*1000/MM;         % Obtain specific gas constant using freestream molar mass [J/kgK]

Freestream.velocity = Velocity - cross(rotrate*[0;0;1],Position); % Freestream velocity compensating for Earth's rotation rate [ECI]- AB


end

