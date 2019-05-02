
close all
clear all
clc

global Vd Z Vm p_amb T_amb  gam a si yi par_IMEP0 par_sp par_fr
global Kthr Kp1 Kp2 Kp3 KT Kfr1 Kfr2 Kfr3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load engine geometric parameters and constant inputs

Vd = 2.4e-3;   % Displacement (m^3)
Z = 4;         % Number of Cylinders
Vm = 5.8e-3;   % Intake Manifold Volume (m^3)
J = 0.0789;    % Mass moment of inertia

p_amb = 1.0121*1e5;
T_amb = 302;
R=288;
gam = 1.35;

P0 = 26431;   % Initial MAP
N0 = 828;     % Initial RPM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters (from steady-state calibration)

a = [1.69e-07,-1.136e-06,6.89e-06];  % Throttle flow model
si = 0.812;   yi = 6.330e3;          % Volumetric efficiency model
par_IMEP0 = [1.2323e-4 2.1256];      % Base IMEP model
par_sp = [-0.0017 -0.0277 1.36];     % Spark timing effects
par_fr = [7.4198e-7 -4.989e-4 11.3]; % Engine friction model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion to Crank Angle Domain

% Coefficients for crank-angle based model
Kthr = p_amb/sqrt(R*T_amb)*sqrt(gam)*sqrt((2/(gam+1))^((gam+1)/(gam-1)));
Kp1 = R*T_amb/Vm;
Kp2 = si*Vd/(4*pi*Vm);
Kp3 = yi*Vd/(4*pi*Vm);
KT = 1e5*Vd/(4*pi);
Kfr1 = (30/pi)^2 * par_fr(1);
Kfr2 = (30/pi) * par_fr(2);
Kfr3 = par_fr(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Equilibrium Condition

setp(1) = 9.81; % Throttle Position  
setp(2) = -25;  % Spark Timing  
setp(3) = 10;   % Load Torque
X0 = fsolve(@(x) linear_Powertrain_Control(x,setp),[26000;24;87]);

P0 = X0(1);           % MAP
T0 = X0(2);           % Tind
N0 = X0(3)*30/pi;     % RPM
disp(sprintf('Equilibrium Value for Pressure: %g Pa',P0));
disp(sprintf('Equilibrium Value for Ind. Torque: %g Nm',T0));
disp(sprintf('Equilibrium Value for Speed: %g r/min',N0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearization

% Coefficients for linearized model as shown in lecture
K1 = Kp1*Kthr*(2*a(1)*setp(1)+a(2));
K2 = Kp1*Kthr*(a(1)*setp(1)^2 + a(2)*setp(1) + a(3));
K3 = Kp2;
Kp = KT*par_IMEP0(1)*(par_sp(1)*setp(2)^2 + par_sp(2)*setp(2) + par_sp(3));    % Pressure multiplier
Kt = KT*(par_IMEP0(1)*X0(1) - par_IMEP0(2)) * (par_sp(1)*setp(2) + par_sp(2)); % Spark Timing multiplier
Kf = 2*Kfr1*X0(3)^2 + Kfr2*X0(3);

tau_D = pi/3;
tau_d = pi/3;
Tsample = pi/3;
d_alpha = 0;
CT1 = 1 + (K3*tau_D)/2;
CT2 = (K3*tau_D)/2 -1;
CT3 = (K2*tau_D)/(2*(X0(3))^2);
CT4 = (K1*tau_D)/(2*X0(3));
CT5 = 1 + (Kf*tau_D)/(2*J*(X0(3))^2);
CT6 = (Kf*tau_D)/(2*J*(X0(3))^2)-1;
CT7 = tau_D/(2*J*X0(3));

% Discretization (solved in Homework 2)
Tau_d = pi/3;   % Discrete sampling time (1 order of engine firing)
Ndel = 3;       % Number of integer delays (in combustion model)
G = tf([CT4*CT7*Kp CT4*CT7*Kp*2 CT4*CT7*Kp],[CT1*CT5 CT1*CT6+CT2*CT5 CT2*CT6 CT3*CT7*Kp CT3*CT7*Kp*2 CT3*CT7*Kp],tau_D,'variable', 'z'  );


M = -(Tau_d)^(0.5).*[CT1 CT3;0 CT5];
N = (Tau_d)^(0.5).*[CT4 0;0 CT7];
O = (Tau_d)^(0.5).*[CT2 CT3;0 CT6];
P = -(Tau_d)^(0.5).*[CT4 0;0 CT7];
Q = P-((O*(inv(M)))*N);
W = inv(M);
Y = W*N;
X = O*(W); 
Phi = [X(1,1) X(1,2) 0 0 Kp*Q(1,2);X(2,1) X(2,2) 0 0 Kp*Q(2,2); W(1,1) W(1,2) 0 0 -Kp*Y(1,2); 0 0 1 0 0; 0 0 0 1 0 ];
Gamma =  [  Q(1,1)         Q(1,2)*Kt;
            Q(2,1)         Q(2,2)*Kt;
            -X(1,1)     -Y(1,2)*Kt;
             0              0;
             0              0;];

H     =  [W(2,1)   W(2,2)     0   0   -Kp*Y(2,2)];       
D     =  [-Y(2,1)   -Y(2,2)*Kt];


Co= ctrb(Phi,Gamma);
uncontrollable_states = length(Phi) - rank(Co);

% %% LQR design
factor = 0.1
Qm = factor*[   0.1           0           0           0           0;
                0           10      0           0           0;
                0           0           0.01           0           0;
                0           0           0           0.01           0;
                0           0           0           0           0.1];
                                
                                
                                
Rm = [  10000  0;
        0   10];

% %% Computation of K 
Ki = [0.005; 0.005];
    
Ksf = dlqr(Phi, Gamma, Qm, Rm); 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Run simulation

[tout,xout,yout] = sim('idle_speed_model_project_LQR');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure('color',[1 1 1],'position',[50 70 700 500])
stairs(omega(:,1), omega(:,2), 'LineWidth', 2);
title({'\fontsize{12} Engine Speed vs Time';'\fontsize{9} Tarey, Aditya'});
ylabel('Engine Speed \omega [r/min]');
xlabel('Time [s]');
grid
% % 
figure('color',[1 1 1],'position',[50 70 700 500])
stairs(al(:,1), al(:,2), 'LineWidth', 2);
title({'\fontsize{12} Throttle Angle vs Time';'\fontsize{9} Tarey, Aditya'});
ylabel('Throttle Angle \alpha [deg]');
xlabel('Time [s]');
grid
% 
figure('color',[1 1 1],'position',[50 70 700 500])
plot(tout,yout(:,2),'linewidth',2)
title({'\fontsize{12} Spark Timing vs Time';'\fontsize{9} Tarey, Aditya'});
ylabel('Spark Timing \theta [bTDC]');
xlabel('Time [s]');
grid
% % 
% % 
% % 
