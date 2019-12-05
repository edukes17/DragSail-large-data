function [ Time_history, State_history ] = Propagation( Windtunnel, State_initial, Orientation_initial, Time_cell, Tolerance, Deorb_Alt, Spacecraft, Surface, outdir, OS)
%% Propagation.m function

%Purpose: Propagate the attitude and trajectory of the spacecraft in orbit.
 
%Created:  Sebastian Tamrazian 12/23/2018
%Modified: Sebastian Tamrazian 2/12/2019 
%Modified: Arly Black 9/24/2019
 
%Acknowledgement: Jay Iuliano for his excellent help and knowledge of orbit
%propagation methods

%Inputs:
%     Windtunnel: 1 or 0 to specify if the mode is active/inactive
%
%     State_initial: Initial user defined state of the spacecraft [1x12]
%                  Initial set of classical orbital elements in the order:
%                     1. Semimajor Axis (m)
%                     2. Eccentricity (-)
%                     3. Inclination (rad)
%                     4. Argument of Perigee (rad)
%                     5. Right Ascension of the Ascending Node (rad)
%                     6. True Anomaly (rad)
%                  Initial Euler angles are in the order:
%                     1. Phi (rad)
%                     2. Theta (rad)
%                     3. Psi (rad)
%                     4. Phidot (rad/s)
%                     5. Thetadot (rad/s)
%                     6. Psidot (rad/s)
%
%     Orientation_initial: [3x1] vector specifying initial euler angle
%     state
%
%     Time_cell: User defined time options for propagation {1x3}
%                     1. Initial date: 'dd-mm-yyyy HH:MM:SS' (string)
%                     2. Propagation time (seconds)
%
%     Tolerance: Specifies solver numerical tolerance 
%
%     Deorb_Alt: User specified deorbit altitude (default 110 km) [m]
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
%     Surface: structure with fields:
%                     1. .temperature : sc surface temperature [K]
%                     2. .reflectance_coef : sc reflectance coefficient - between 0 and 1 []
%                     3. .specular_coef : sc specular reflectance coefficient - between 0 and 1 []
%                     4. .front_emiss_coef : front emissivity coefficient - between 0 and 1 []
%                     5. .back_emiss_coef : back emissivity coefficient - between 0 and 1 []
%                     6. .front_nonLamb_coef : front non-Lambertian coefficient - between 0 and 1 []
%                     7. .back_nonLamb_coef : back non-Lambertian coefficient - between 0 and 1 []
%
%     outdir: specifies output file path for simulation results
%
%     OS: specifies operating system type: PC or Mac (impacts file directory storage syntax) 


% Supporting Functions: 
%    Perturbations_Aerodynamics.m
%    Perturbations_Gravitational.m
%    Perturbations_SolarRadiationPressure.m
%    AtmosphericProperties.m
%    eclipseCheck.m
%    SunPosition.m
%    Convert_BODY2ECI.m
%    Convert_COES2ECI.m
%    Convert_EA2EP.m
%    PlotResults.m
                                  
%% Initialize Variables

% Earth Properties
EARTH.SGP                   = 3.9860044189e14;          % Std. grav. param. [m^3/s^2]
EARTH.MASS                  = 5.97e24;                  % Mass [kg]
EARTH.EQRADIUS              = 6378.1363e3;              % Equatorial radius [m]
EARTH.PORADIUS              = 6051.8e3;                 % Polar radius [m]
EARTH.J2                    = 1082.63e-6;               % J2 oblateness parameter [-]


% Initial orbital elements (classical)
a_initial                   = State_initial(1);         % Semi-major axis [m] 
e_initial                   = State_initial(2);         % Eccentricity []
i_initial                   = State_initial(3);         % Inclination [rad]
w_initial                   = State_initial(4);         % Arguement of periapsis [rad]
o_initial                   = State_initial(5);         % Raan [rad]
t_initial                   = State_initial(6);         % True anomaly [rad]

% Initial attitude state
Angularrates_initial        = State_initial(7:9);       % [rad/s]
Quaternions_initial         = State_initial(10:13);     % [q1 q2 q3 q4]

% Initial radius and velocity vectors from classical orbital elements
AngularMomentum_initial = sqrt(a_initial*(1-e_initial^2)*EARTH.SGP);    % angular momentum [m2/s]
[Position_initial,Velocity_initial] = Convert_COES2ECI(AngularMomentum_initial,t_initial,w_initial,o_initial,i_initial,e_initial,EARTH.SGP); 


% Initial state vectors:
ThreeDOF_initial            = [Quaternions_initial,Angularrates_initial]; % Windtunnel mode
SixDOF_initial              = [Position_initial', Velocity_initial', Quaternions_initial, Angularrates_initial]; % Orbital mode

% Time Options
Time_epoch                  = Time_cell{1};    % 'mm-dd-yyyy HH:MM:SS'
Time_propagation            = Time_cell{2};    % [seconds]

%% Propagate 3DOF or 6DOF depending on if Windunnel mode is active

if ~Windtunnel 
    %Runs if Orbital mode (6DOF) is enabled
    
    % Set options for specified tolerance and check for deorbit
    options  = odeset('Abstol', Tolerance, 'Reltol', Tolerance, 'Events', @Deorbit);
    
    % Find position of sun
    [R_sun] = SunPosition(Time_epoch);
 
    %define timeseries to save parameters    
    global ts1
    global ts2
    global ts3
    global ts4
    global ts5
    global ts6
    
    ts1 = timeseries('AeroTorques');
    ts2 = timeseries('AeroForces');
    ts3 = timeseries('GravTorques');
    ts4 = timeseries('SolarTorques');
    ts5 = timeseries('SolarForces');
    ts6 = timeseries('SunBody');
    
    
    %Propagate
    [Time_history,State_history] = ode113(@(Time, State) Orbitalcase(Time, State, Spacecraft, Surface, EARTH, Time_epoch, R_sun),...
                                            [0 Time_propagation], SixDOF_initial, options);

    State_history = real(State_history);
    
    % Pre-allocate variable size
    zeros3 = zeros(3,length(State_history));
    Torques_aero = zeros3;
    Forces_aero  = zeros3;
    Torques_srp  = zeros3;
    Forces_srp   = zeros3;
    Torques_grav = zeros3;
    Forces_grav  = zeros3;
    alpha        = zeros3(1,:);
    beta         = zeros3(1,:);
    aoa         = zeros3(1,:);
    Vrel_body    = zeros3;
    SB           = zeros3;
    Cdrag        = zeros3(1,:);
    Clift        = zeros3(1,:);
    
    % Calculate forces and torques for given state history
    for nn = 1:(length(State_history))
        r = State_history(nn, 1:3);             % position [m]  
        v = State_history(nn, 4:6);             % velocity [m/s]
        q = State_history(nn, 7:10);            % quaternions, q4 scalar [none]
        
               
        Freestream = AtmosphericProperties(EARTH, Time_epoch, State_history(nn,1:3)', State_history(nn,4:6)'); % freestream properties
        Euler_angles = Convert_EP2EA(q(1:4));                   % convert from quaternions to euler angles
        BODY2ECI = Convert_BODY2ECI(Euler_angles);              % convert from body to ECI frame
        Shadow = eclipseCheck(R_sun,r',EARTH);                  % check if s/c is in eclipse (for SRP calcs)
        vr_body = BODY2ECI'*Freestream.velocity;
        alpha(1,nn) = atan2(vr_body(3),vr_body(1));             % calculate angle of attack [rad]
        beta(1,nn) = atan2(vr_body(2),vr_body(1));              % calculate side slip angle [rad]
        vrhat = vr_body/norm(vr_body);
        aoa(1,nn) = acos(dot(vrhat,[1,0,0]));              % total angle of attack [rad] (angle bw vr and x axis)
        
        % Calculate perturbations
        [~,Torque_aero,Force_aero,Cd,Cl] = Perturbations_Aerodynamics(Spacecraft, Freestream, Surface, BODY2ECI); 
        [~, Torque_grav,Force_asph] = Perturbations_Gravitational(r', EARTH, Spacecraft, BODY2ECI);
        [~, Torque_srp,Force_srp,R_sun_hat] = Perturbations_SolarRadiationPressure(Spacecraft, r', R_sun, Surface, Shadow, BODY2ECI);
        
        % Store values
%        Torques_aero(:,nn) = Torque_aero';
%        Forces_aero(:,nn)  = Force_aero';
%        Torques_grav(:,nn) = Torque_grav';
        Forces_grav(:,nn)  = Force_asph';
%        Torques_srp(:,nn)  = Torque_srp';
%        Forces_srp(:,nn)   = Force_srp';
        Vrel_body(:,nn)    = vr_body;
%        SB(:,nn)           = R_sun_hat';          % Sun-in-body vector
        Cdrag(nn,1)        = Cd;
        Clift(nn,1)        = Cl;
        
        
    end
       
    Euler_history = rad2deg(Convert_EP2EA(State_history(:,7:10)));
    ang_of_attack = rad2deg(alpha');
    sideslip = rad2deg(beta');
    aoa_tot = rad2deg(aoa');
    
    vr = Vrel_body';
    F_grav = Forces_grav';
    T_grav = Torques_grav';
    F_srp = Forces_srp';
    T_srp = Torques_srp';
    Sun_body = SB';
    F_aero = Forces_aero';
    T_aero = Torques_aero';

    
    %% Output results
    
    T1 = table(Time_history,State_history,ang_of_attack,sideslip,aoa_tot,vr);
    T2 = table(Time_history,T_grav,F_grav,Euler_history);
    T3 = table(Time_history,T_srp,F_srp,Euler_history,Sun_body);
    T4 = table(Time_history,T_aero,F_aero,Cdrag,Clift);
    
    % Check if user OS is PC or Mac, then output results to indicated folder
    if ispc 
        % For PC
        writetable(T1,[pwd '\' outdir '\StateHistory.csv']);
        writetable(T2,[pwd '\' outdir '\grav.csv']);
        writetable(T3,[pwd '\' outdir '\SRP.csv']);
        writetable(T4,[pwd '\' outdir '\aero.csv']);

    elseif ismac
        % For Mac
        writetable(T1,[pwd '/' outdir '/StateHistory.csv']);
        writetable(T2,[pwd '/' outdir '/grav.csv']);
        writetable(T3,[pwd '/' outdir '/SRP.csv']);
        writetable(T4,[pwd '/' outdir '/aero.csv']);
    end
    
    % Plot results
    PlotResults(Time_history, State_history, Euler_history, Orientation_initial, Torques_aero, Forces_aero, Torques_grav, Forces_grav,...
        Torques_srp, Forces_srp, alpha, beta, aoa, EARTH, R_sun, Sun_body, outdir);
   
    figure(12)
    plot(ts1);
    legend('Tx','Ty','Tz');
    ylabel('Torque [Nm]');xlabel('Time [sec]');
    
    T5 = table(ts1.time,ts1.data,'VariableNames',{'Time','TA'});
    writetable(T5,[pwd '/' outdir '/aeroTorqueForce.csv']);

end
end


%% Orbital Formulation ODE function

 function [Statederivative] = Orbitalcase(Time, State, Spacecraft, Surface, EARTH, Time_epoch, R_sun)
          
   global ts1
   global ts2
   global ts3
   global ts4
   global ts5
   global ts6
   
 % Translational state
    r = State(1:3); %position [m]  
    v = State(4:6); %velocity [m/s]  
 
    % Rotational state
    q = real(State(7:10));                          % quaternions, q(4) is the scalar component
    w = State(11:13);                               % angular rates [rad/s]
    Euler_angles = Convert_EP2EA(q(1:4)');
    [BODY2ECI] = Convert_BODY2ECI(Euler_angles);    % body to ECI coordinate conversion
        
    % Obtain freestream conditions and earth shadow for current state
    [Freestream] = AtmosphericProperties(EARTH, Time_epoch, r, v);
    [Shadow] = eclipseCheck(R_sun,r,EARTH);
    
    % Find  perturbations
    [Accel_aerodynamic,Torque_aerodynamic,Force_aerodynamic, Cd, Cl] = Perturbations_Aerodynamics(Spacecraft, Freestream, Surface, BODY2ECI);
    %Save the results
    Torque_aero = transpose(Torque_aerodynamic);
    Force_aero = transpose(Force_aerodynamic);
    ts1 = addsample(ts1,'Data',Torque_aero,'Time',Time);
    ts2 = addsample(ts2,'Data',Force_aero,'Time',Time);
    
    [Accel_aspherical, Torque_gravity,~] = Perturbations_Gravitational(r, EARTH, Spacecraft, BODY2ECI);
    Torque_grav = transpose(Torque_gravity);
    ts3 = addsample(ts3,'Data',Torque_grav,'Time',Time);
    
    [Accel_srp, Torque_srp,Force_srp,R_sun_hat] = Perturbations_SolarRadiationPressure(Spacecraft, r, R_sun, Surface, Shadow, BODY2ECI);
    Torque_srp = transpose(Torque_srp);
    Force_srp = transpose(Force_srp);
    Sun_body = transpose(R_sun_hat);
    ts4 = addsample(ts4,'Data',Torque_srp,'Time',Time);
    ts5 = addsample(ts5,'Data',Force_srp,'Time', Time);
    ts6 = addsample(ts6,'Data',Sun_body,'Time', Time);
        
    % Sum Perturbations
    Accelerations = Accel_aerodynamic + Accel_aspherical + Accel_srp;
    Torques = Torque_aerodynamic + Torque_gravity + Torque_srp;
    
  % Translation :
    a = -EARTH.SGP/norm(r)^3*r + Accelerations;
       
  % Rotation :
    % Dynamic EOMs:
    % Torque = I*wdot + (w x H)
    H = Spacecraft.inertia*w;           % angular momentum = inertia tensor * angular velocity [kg*m2/s]
    BKE = cross(w,H);                   % w x H [kg*m2/s2 = Nm]
    I_wdot = Torques-BKE;               % I*wdot = Torque - (w x H) [Nm]
    wd = I_wdot./Spacecraft.inertia;    % wdot matrix
    wdot = [wd(1,1);wd(2,2);wd(3,3)];   % angular rate derivatives [rad/s2]

    % Kinematic EOMs (quaternion formulation):
    qdot(1,1) = (w(1)*q(4)-w(2)*q(3)+w(3)*q(2))/2;
    qdot(2,1) = (w(1)*q(3)+w(2)*q(4)-w(3)*q(1))/2;
    qdot(3,1) =(-w(1)*q(2)+w(2)*q(1)+w(3)*q(4))/2;
    qdot(4,1) =-(w(1)*q(1)+w(2)*q(2)+w(3)*q(3))/2;

    % 6DOF state derivative
    Statederivative = [v;a;qdot;wdot];
 
    % Print Status
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b') % delete previous time
    fprintf('%13.2f',Time) % print current time

    
 end

function [value,isterminal,direction] = Deorbit(~, State, ~, ~, ~, ~, ~, ~)    

    r = State(1:3); % position [m]
    global Deorb_Alt
    
    %determine altitude
    Altitude = norm(r)-6378100;
    if Altitude <= Deorb_Alt
        check = 1;
    else
        check = 0;
    end
    
    %if Altitude is less than "Deorb_Alt" m, propagation stops
    value = 1 - check;
    isterminal = 1;
    direction = 0;
end
 
