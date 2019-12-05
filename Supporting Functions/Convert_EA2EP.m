function [Quaternions] = Convert_EA2EP(Euler_angles)
%% Convert_EA2EP Function

%Purpose: Convert euler angles to euler parameters (quaternions). This
%formulation is specific to a body 3-2-1 rotation sequence.
 
%Created:  Sebastian Tamrazian 1/7/2019
 
%Inputs:    Euler_angles [Mx3] matrix [yaw pitch roll], [rad]
%Outputs:   Quaternions [Mx4] matrix [q1 q2 q3 q4], with q4 as the scalar


%Scalar component        
Quaternions(:,4) =  sin(Euler_angles(:,1)/2).*sin(Euler_angles(:,2)/2).*sin(Euler_angles(:,3)/2)...
                   +cos(Euler_angles(:,1)/2).*cos(Euler_angles(:,2)/2).*cos(Euler_angles(:,3)/2);

%Vector components
Quaternions(:,1) =  sin(Euler_angles(:,1)/2).*cos(Euler_angles(:,2)/2).*cos(Euler_angles(:,3)/2)...
                   -cos(Euler_angles(:,1)/2).*sin(Euler_angles(:,2)/2).*sin(Euler_angles(:,3)/2);
               
Quaternions(:,2) =  cos(Euler_angles(:,1)/2).*sin(Euler_angles(:,2)/2).*cos(Euler_angles(:,3)/2)...
                   +sin(Euler_angles(:,1)/2).*cos(Euler_angles(:,2)/2).*sin(Euler_angles(:,3)/2);

Quaternions(:,3) =  cos(Euler_angles(:,1)/2).*cos(Euler_angles(:,2)/2).*sin(Euler_angles(:,3)/2)...
                   -sin(Euler_angles(:,1)/2).*sin(Euler_angles(:,2)/2).*cos(Euler_angles(:,3)/2);
               

% %Scalar component        
% Quaternions(:,4) =  sin(Euler_angles(:,3)/2).*sin(Euler_angles(:,2)/2).*sin(Euler_angles(:,1)/2)...
%                    +cos(Euler_angles(:,3)/2).*cos(Euler_angles(:,2)/2).*cos(Euler_angles(:,1)/2);
% 
% %Vector components
% Quaternions(:,1) =  sin(Euler_angles(:,3)/2).*cos(Euler_angles(:,2)/2).*cos(Euler_angles(:,1)/2)...
%                    -cos(Euler_angles(:,3)/2).*sin(Euler_angles(:,2)/2).*sin(Euler_angles(:,1)/2);
%                
% Quaternions(:,2) =  cos(Euler_angles(:,3)/2).*sin(Euler_angles(:,2)/2).*cos(Euler_angles(:,1)/2)...
%                    +sin(Euler_angles(:,3)/2).*cos(Euler_angles(:,2)/2).*sin(Euler_angles(:,1)/2);
% 
% Quaternions(:,3) =  cos(Euler_angles(:,3)/2).*cos(Euler_angles(:,2)/2).*sin(Euler_angles(:,1)/2)...
%                    -sin(Euler_angles(:,3)/2).*sin(Euler_angles(:,2)/2).*cos(Euler_angles(:,1)/2);
               
             
end

