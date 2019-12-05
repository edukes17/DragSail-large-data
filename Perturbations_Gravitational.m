function [Accel_Aspherical, Torque_Gravity, Force_Gravity] = Perturbations_Gravitational(Position, EARTH, Spacecraft, BODY2ECI)
%% Perturbations_Gravitational.m function
%Purpose: Determine gravitational forces and torques on the spacecraft for
%given geometry and state

% Created:  Sebastian Tamrazian 2/1/2019
% Modified: Arly Black 8/15/2019

% Inputs:

%     Position: spacecraft position vector in ECI [m]
%
%     EARTH: structure with fields:
%                     1. .SGP : Earth standard gravitational parameter [m^3/s^2]
%                     2. .MASS : Earth mass [kg]
%                     3. .EQRADIUS : Earth equatorial radius [m]
%                     4. .PORADIUS : Earth polar radius [m]
%                     5. .J2 : first zonal harmonic coefficient of Earth [-]
%                     6. .rotrate : Earth rotation rate [rad/s]
%
%     Spacecraft: structure with fields:
%                     1. .mass : sc mass [kg]
%                     2. .inertia: spacecraft inertia [kg*m^2]
%                     3. .cm : spacecraft cm relative to geom origin [m]
%                     4. .centroids : centroid locations of sc surface elements [m]
%                     5. .normals : normals of sc surface elements (unit vectors) []
%                     6. .areas : areas of sc surface elements [m^2]
%                     7. .tangents : tangents of sc surface elements (unit vectors)[]
%
%
%     ROT2ECI: [3x3] Rotation matrix in form [ECI] = ROT2ECI[ROT]
%
%     BODY2ROT: [3x3] Rotation matrix in form [ROT] = BODY2ROT[BODY]

%% Orbital perturbation:

% First zonal harmonics term, J2, is accounted for with a simplified
% acceleration model. Other zonal harmonics terms, as well as sectoral and
% tessoral harmonics terms, are omitted in this simple formulation.

Distance = norm(Position);
    
Coef_J2 = -3/2*EARTH.J2*EARTH.SGP*EARTH.EQRADIUS^2/Distance^5;
Accel_J2(1,1) = Coef_J2*(1-5*(Position(3)/Distance)^2)*Position(1);
Accel_J2(2,1) = Coef_J2*(1-5*(Position(3)/Distance)^2)*Position(2);
Accel_J2(3,1) = Coef_J2*(3-5*(Position(3)/Distance)^2)*Position(3);

Accel_Aspherical = Accel_J2; % Accelerations [m/s2]
Force_Gravity = Accel_Aspherical*Spacecraft.mass; % [N = kg*m/s2]

%% Attitude perturbation:

% The gravity gradient torque is found using the spacecraft inertia
% properties. 

Position_Body = BODY2ECI'*Position;
   
I11 = Spacecraft.inertia(1,1);
I12 = Spacecraft.inertia(1,2);
I13 = Spacecraft.inertia(1,3);      % In = [I11 I12 I13
I21 = Spacecraft.inertia(2,1);      %       I21 I22 I23
I22 = Spacecraft.inertia(2,2);      %       I31 I32 I33]
I23 = Spacecraft.inertia(2,3);
I31 = Spacecraft.inertia(3,1);
I32 = Spacecraft.inertia(3,2);
I33 = Spacecraft.inertia(3,3);


Nad = -Position_Body./norm(Position_Body); % Nadir unit vector expressed in body-fixed frame

NIN = [Nad(2)*Nad(3)*(I13+I23+I33-I12-I22-I32);     % = Nad x (In)(Nad)
       Nad(3)*Nad(1)*(I11+I21+I31-I13-I23-I33);
       Nad(1)*Nad(2)*(I12+I22+I32-I11-I21-I31)];

Torque_Gravity = 3*EARTH.SGP/(norm(Position_Body))^3*NIN; % Gravity gradient torque vector in body frame [Nm]

end