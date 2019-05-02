function y = linear_ECE753(x,par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Example of function to be used with "fsolve".
%
% This function implements a set of nonlinear algebraic equations that
% represent the steady-state equilibrium condition for the linearized
% system. Used with "fsolve", this function allows for determining the
% values of manifold pressure and engine speed at equilibrium, for a given
% combination of throttle opening and disturbance torque.
%
%    * INPUTS: x (equilibrium values of manifold pressure and speed)
%              par (equilibrium values of throttle and dist. torque)
%    * OUTPUT: y (2 nonlinear equations obtained during linearization)
%
%
% (C) Marcello Canova, 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define global variables

global a Kthr Kp1 Kp2 Kp3 KT Kfr1 Kfr2 Kfr3 par_IMEP0 par_sp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables:
% x = [p0 T0 w0]
% par = [alpha, spark, Tload]

% Write Model Equations in SS Conditions (to be imposed equal to zero)
y = [Kp1*Kthr*(a(1)*par(1)^2 + a(2)*par(1) +a(3)) - Kp2*x(1)*x(3) + Kp3*x(3);
     x(2) - KT*(par_IMEP0(1)*x(1) - par_IMEP0(2)) * (par_sp(1)*par(2)^2 + par_sp(2)*par(2) + par_sp(3));
     x(2) - par(3) - Kfr1*x(3)^2 - Kfr2*x(3) - Kfr3];
 
