function [Euler_angles] = Convert_EP2EA(Quaternions)
%% Convert_EP2EA Function

%Purpose: Convert euler parameters (quaternions) to euler angles
 
%Created:  Sebastian Tamrazian 1/7/2019

%Inputs:   Quaternions [1x4] vector [q1 q2 q3 q4], with q4 as the scalar
%component
%Outputs:  Euler_angles [1x3] vector [yaw pitch roll], [rad]

%Variables: 
%   phi  : yaw [rad]
%   theta: pitch [rad]
%   psi  : roll [rad]
%   tol  : tolerance for calculation in fixing gimbal lock


q1 = Quaternions(:,1);
q2 = Quaternions(:,2);
q3 = Quaternions(:,3);
q4 = Quaternions(:,4);
tol = .0001;

e1 = 2*(q4.*q1 + q2.*q3);
e2 = 1 - 2*(q1.^2 + q2.^2);
e3 = 2*(q4.*q3 + q1.*q2);
e4 = 1 - 2*(q2.^2 + q3.^2);

% Euler angles calculation
phi = atan2(2*(q4.*q1 + q2.*q3),1 - 2*(q1.^2 + q2.^2));
theta = asin(2*(q4.*q2 - q3.*q1));
psi = atan2(2*(q4.*q3 + q1.*q2),1 - 2*(q2.^2 + q3.^2));

% Gimbal lock (singularities) adjustment
phi(theta > pi/2-tol) = 0;
theta(theta > pi/2-tol) = pi/2;
psi(theta > pi/2-tol) = 2*atan2(q1(theta > pi/2-tol),q4(theta > pi/2-tol));

phi(theta < -pi/2+tol) = 0;
theta(theta < -pi/2+tol) = -pi/2;
psi(theta < -pi/2+tol) = -2*atan2(q1(theta < -pi/2+tol),q4(theta < -pi/2+tol));

%compile results
Euler_angles = [phi,theta,psi]; %[yaw,pitch,roll]
end

