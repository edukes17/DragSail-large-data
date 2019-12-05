function [Shadow] = eclipseCheck(R_sun, r, EARTH)

% Purpose: Determine whether spacecraft is in sunlight or eclipse - affects
% solar radiation pressure forces and torques

% Incorporated into AeroCode  ED 10/14/19

% Inputs:
% r                 s/c position vector in ECI frame (m)
% R_sun             Sun position vector in ECI frame (m)
% EARTH.EQRADIUS    Earth radius (m)

% Output:
% Shadow            Binary value describing whether s/c in in sunlight (1)
%                   or eclipse(0)


% Re = 6378.1363e3;   %Radius of Earth (m)
Re = EARTH.EQRADIUS;   %Radius of Earth (m)

R_sunHat = R_sun/norm(R_sun); %Sun unit vector in ECI frame
rHat = r/norm(r); %Earth unit vector in ECI frame
uHat = cross(R_sunHat,rHat); %rSunHat x rHat
pHat = cross(uHat,R_sunHat); %uHat x rSunHat, perpendicular to Earth-Sun line in Sun-Earth-s/c plane

SEPangle = acos(dot(R_sunHat,rHat)); %Sun-Earth-s/c angle
dot_r_pHat = dot(r,pHat); %dot product between r and pHat; projection of r onto pHat direction

if SEPangle <= pi/2 || dot_r_pHat > Re
    Shadow = 1; %not in eclipse
else
    Shadow = 0; % In eclipse
end

end