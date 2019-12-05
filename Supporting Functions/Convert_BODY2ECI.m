function [BODY2ECI] = Convert_BODY2ECI(Euler_angles)
%% Convert_BODY2ROT Function

%Purpose: Convert from spacecraft body frame to Earth
%centered inertial coordinates (ECI or IJK). Uses a 3 2 1 (yaw pitch roll)
%rotation sequence

%Created:  Arly Black 9/20/2019
 
%Inputs:    Euler angles in the order roll, pitch, yaw

%Outputs:   Rotation matrix:
%                   [x;           [x1 x2 x3;     [X_r;
%  ECI2BODY          y;     =      y1 y2 y3;  *   Y_r;
%                    z]            z1 z2 z3]      Z_r]

%                  [X_r;          [x1 y1 z1;      [x;
%  BODY2ECI         Y_r;    =      x2 y2 z2;  *    y;
%                   Z_r]           x3 y3 z3]       z]

phi = Euler_angles(1);      % roll
theta = Euler_angles(2);    % pitch
psi = Euler_angles(3);      % yaw
       
sph = sin(phi);
cph = cos(phi);
st = sin(theta);
ct = cos(theta);
sp = sin(psi);
cp = cos(psi);

ECI2BODY = [ct*cp, ct*sp, -st;...
            sph*st*cp-cph*sp, sph*st*sp+cph*cp, sph*ct;...
            cph*st*cp+sph*sp, cph*st*sp-sph*cp, cph*ct];      % Direction cosine matrix
        
BODY2ECI = ECI2BODY';

end

