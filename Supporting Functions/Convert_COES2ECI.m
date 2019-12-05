%{
Calulate ECI position and velocity from classical orbital elements

Created: Jay Iuliano 9/19/2017
Modified: Jay Iuliano 4/4/2018 - Modified to accept vectors

Inputs:
h            -> [nx1] Orbital angular momentum (m^2/s)
theta        -> [nx1] True anomaly (rad)
w            -> [nx1] Argument of perigee (rad)
raan         -> [nx1] Right ascension of ascending node (rad)
inc          -> [nx1] Inclination (rad)
ecc          -> [nx1] Eccentricity (none)
MU           -> [nx1] Gravitational constant of primary body (m^3/s^2)

Outputs:
r            -> [3xn] Radius vector in ECI frame (m)
v            -> [3xn] Velocity vector in ECI frame (m/s)
%}

function [r,v] = Convert_COES2ECI(h,theta,w,raan,inc,ecc,MU)

% Precalculate
ct = cos(theta); st = sin(theta);
cw = cos(w); sw = sin(w);
cr = cos(raan); sr = sin(raan);
ci = cos(inc); si = sin(inc);

if numel(theta) == 1
    A = h*h/MU/(1 + ecc*ct);
    B = MU/h;

    % Perifocal Coordinates
    x_peri = A*ct; 
    y_peri = A*st; 
    % z_peri = zeros(length(h),1);

    vx_peri = -B*st; 
    vy_peri = B*(ecc + ct); 
    % vz_peri = zeros(length(h),1);

    % Preallocate DCM values for speed
    DCM_peri2ECI_11 = (cw*cr - sw*ci*sr);
    DCM_peri2ECI_12 = (-sw*cr - cw*ci*sr);

    DCM_peri2ECI_21 = (cw*sr + sw*ci*cr);
    DCM_peri2ECI_22 = (-sw*sr + cw*ci*cr);

    DCM_peri2ECI_31 = si*sw;
    DCM_peri2ECI_32 = si*cw;

    % Inertial position and velocity
    x = DCM_peri2ECI_11*x_peri + DCM_peri2ECI_12*y_peri;
    y = DCM_peri2ECI_21*x_peri + DCM_peri2ECI_22*y_peri;
    z = DCM_peri2ECI_31*x_peri + DCM_peri2ECI_32*y_peri;
    r = [x;y;z];

    vx = DCM_peri2ECI_11*vx_peri + DCM_peri2ECI_12*vy_peri;
    vy = DCM_peri2ECI_21*vx_peri + DCM_peri2ECI_22*vy_peri;
    vz = DCM_peri2ECI_31*vx_peri + DCM_peri2ECI_32*vy_peri;
    v = [vx;vy;vz];
else
    A = h.*h./MU./(1 + ecc.*ct);
    B = MU./h;

    % Perifocal Coordinates
    x_peri = A.*ct; 
    y_peri = A.*st; 
    % z_peri = zeros(length(h),1);

    vx_peri = -B.*st; 
    vy_peri = B.*(ecc + ct); 
    % vz_peri = zeros(length(h),1);

    % Preallocate DCM values for speed
    DCM_peri2ECI_11 = (cw.*cr - sw.*ci.*sr);
    DCM_peri2ECI_12 = (-sw.*cr - cw.*ci.*sr);

    DCM_peri2ECI_21 = (cw.*sr + sw.*ci.*cr);
    DCM_peri2ECI_22 = (-sw.*sr + cw.*ci.*cr);

    DCM_peri2ECI_31 = si.*sw;
    DCM_peri2ECI_32 = si.*cw;

    % Inertial position and velocity
    x = DCM_peri2ECI_11.*x_peri + DCM_peri2ECI_12.*y_peri;
    y = DCM_peri2ECI_21.*x_peri + DCM_peri2ECI_22.*y_peri;
    z = DCM_peri2ECI_31.*x_peri + DCM_peri2ECI_32.*y_peri;
    r = [x,y,z]'; 

    vx = DCM_peri2ECI_11.*vx_peri + DCM_peri2ECI_12.*vy_peri;
    vy = DCM_peri2ECI_21.*vx_peri + DCM_peri2ECI_22.*vy_peri;
    vz = DCM_peri2ECI_31.*vx_peri + DCM_peri2ECI_32.*vy_peri;
    v = [vx,vy,vz]';
end