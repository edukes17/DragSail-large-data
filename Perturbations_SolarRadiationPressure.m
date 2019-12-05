function [ Acceleration_srp, Torque_srp, Force_srp, R_sun_hat ] = Perturbations_SolarRadiationPressure( Spacecraft, Position, R_sun, Surface, EarthShadow, BODY2ECI)
%% SolarRadiationPressure.m function

% Purpose: Determine SRP forces and torques on the spacecraft for given geometry and state

% Created:  Sebastian Tamrazian 2/6/2019
% Updated:  Eileen Dukes 9/23/19 to incorporate Solar Sail Force Model from McInnes, "Solar Sailing"
% Updated:  Eileen Dukes 11/25/19 to accommodate clear sail

% Inputs:
%     Spacecraft: structure with fields:
%                     1. .mass : sc mass [kg]
%                     2. .inertia: spacecraft inertia [kg*m^2]
%                     3. .cm : spacecraft cm relative to geom origin [m]
%                     4. .centroids : centroid locations of sc surface elements [m]
%                     5. .normals : normals of sc surface elements (unit vectors) []
%                     6. .areas : areas of sc surface elements [m^2]
%                     7. .tangents : tangents of sc surface elements (unit vectors)[]
%
%     Position: spacecraft position vector in ECI [m]
%
%     R_sun: Sun position vector in ECI [m]
%
%     Surface: structure with fields:
%                     1. .temperature : sc surface temperature [K]
%                     2. .reflectance_coef : sc reflectance coefficient [] non-dimensional coefficient between 0 and 1 where 1=total reflector
%                     3. .specular_coef = Specular Reflectance  - non-dimensional coefficient between 0 and 1
%                     4. .front_emiss_coef =  front emissivity - non-dimensional coefficient between 0 and 1
%                     5. .back_emiss_coef =  Back emissivity - non-dimensional coefficient between 0 and 1
%                     6. .front_nonLamb_coef =  Front non-Lambertian coefficient - non-dimensional coefficient between 0 and 1
%                     7. .back_nonLamb_coef =  Back non-Lambertian coefficient - non-dimensional coefficient between 0 and 1
%                     8. .transmissivity_coef = Transmissivity for use with clear sail
%                     9. .Refl                = Specify clear or reflecting surface 1=Refl, 0=Clear
%             *This implementation assumes that all of the surfaces have the same coefficients These need to be included in the inputs
%

%     EarthShadow: 0 = spacecraft is in eclipse (SRP perturbations = 0)
%                  1 = full sunlight
%
%     ROT2ECI: [3x3] Rotation matrix in form [ECI] = ROT2ECI[ROT]
%
%     BODY2ROT: [3x3] Rotation matrix in form [ROT] = BODY2ROT[BODY]

%  Outputs:
%   Solar force vector in the BODY frame
%   Solar torque vector in the BODY frame
%   Solar acceleration vector in the ECI frame

%If the spacecraft is in eclipse, then the solar forces are zero
if EarthShadow == 0
    Acceleration_srp = [0;0;0];
    Torque_srp = [0;0;0];
    Force_srp = [0;0;0];
    R_sun_hat = [0;0;0];
    return
end


p_srp = 4.571176e-6; %pressure from sun photons at Earth orbit [N/m^2] 

% introduce new coefficients for clear sail
ab = 1-Surface.transmissivity_coef;
p_spec = Surface.reflectance_coef * Surface.specular_coef;   %specular reflection coefficient (rs)
d = Surface.reflectance_coef * (1-Surface.specular_coef);   %diffuse reflections coefficient r(1-s)
Refl = Surface.Refl;


% Rotate sun vector into the body frame to determine incidence angle with
% panels
R_sat_sun = R_sun-Position;                 % Vector from satellite to sun [ECI]
R_sun_body = BODY2ECI'*R_sat_sun;
R_sun_hat = R_sun_body/norm(R_sun_body);    % Unit vector of sun vector [BODY] 
Sun_in = -R_sun_hat;                        % Incoming sun vector to match McInnes Fig 2.7


Force_srp_vector = zeros(size(Spacecraft.normals))';
Torque_srp_vector = zeros(size(Spacecraft.normals))';

% Calculate the emissivity and Non-Lambertian term which is only dependent on surface properties 
Lamb_emiss = ((Surface.front_emiss_coef*Surface.front_nonLamb_coef)-(Surface.back_emiss_coef*Surface.back_nonLamb_coef))...
            /(Surface.front_emiss_coef+Surface.back_emiss_coef);
 
% Iterate through all surface elements
for jj = 1:length(Spacecraft.areas)

    % Obtain element information - need to add calculation of the element
    % tangent vector in Spacecraft ED 9/24/19
    Normal_jj = Spacecraft.normals(jj,:);
    Area_jj = Spacecraft.areas(jj);
    Centroid_jj = Spacecraft.centroids(jj,:);
    Tangent_jj = Spacecraft.tangents(jj,:);
    
    % Determine angle offset from sun vector- fixed error to use acos
    % deleted negative sign, Sun vector already in body frame - ED 9/24/19
    Alpha = acos(dot(Sun_in,Normal_jj));
    
    % Simple shadowing check - If the surface is pointed away from the Sun consider it shadowed ED-9/24/19
    % Check to see if the surface normal is greater than 89.9 deg from the Sun vector 
    test= abs(Alpha);
        if test>deg2rad(89.9)
              
           PanelShadow=0; % shadowed
        else
           PanelShadow=1; % not shadowed
        end
    
    % Find surface cell force and moment
    % equations depend on whether the material has a reflective coating
    
    if Refl == 1   %reflective coating
    
    % Calculate the force normal to the surface element
    % Remove Negative sign, normal vector switched
        Force_normal_srp_jj = p_srp*Area_jj*((1+Surface.reflectance_coef*Surface.specular_coef)...
             *(cos(Alpha)*cos(Alpha))+Surface.front_nonLamb_coef*(1-Surface.specular_coef)*Surface.reflectance_coef...
             *cos(Alpha)+(1-Surface.reflectance_coef)*cos(Alpha)*Lamb_emiss)*Normal_jj*PanelShadow*EarthShadow; %corrected last term
   
    % Calculate the force tangental to the surface element
         Force_tang_srp_jj = (p_srp*Area_jj*(1-Surface.reflectance_coef*Surface.specular_coef)*cos(Alpha)*sin(Alpha))*...
         Tangent_jj*PanelShadow*EarthShadow;
    else
        %Calculate the force normal to the surface element for a clear sail
         Force_normal_srp_jj = p_srp*Area_jj*(((ab+p_spec)*cos(Alpha)*cos(Alpha))+(Surface.front_nonLamb_coef*d*cos(Alpha))...
             +(1-Surface.reflectance_coef)*cos(Alpha)*Lamb_emiss)*...
             Normal_jj*PanelShadow*EarthShadow;
         %Calculate the force tangent for a clear sail
         Force_tang_srp_jj = p_srp*Area_jj*((ab-p_spec)*cos(Alpha)*sin(Alpha))*Tangent_jj*PanelShadow*EarthShadow;
    end
    % Sum the forces
    Force_srp_jj = Force_normal_srp_jj+Force_tang_srp_jj;
     
    %Calculate the torque 
    Torque_srp_jj = cross(Centroid_jj, Force_srp_jj);
    
    %Add force and moment to total
    Force_srp_vector(:,jj) = Force_srp_jj;
    Torque_srp_vector(:,jj) = Torque_srp_jj;

end

Force_srp = sum(Force_srp_vector,2); % [N] - Body frame
Torque_srp = sum(Torque_srp_vector,2); % [Nm] - Body frame
Acceleration_srp = (BODY2ECI*Force_srp)/Spacecraft.mass; % [m/s2] - ECI frame
end

