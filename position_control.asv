% 
% Copyright (c) 2020 Flight Dynamics and Control Lab
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the 
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
%  in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY,FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%The code is modified according to Control of Complex Maneuvers for a Quadrotor UAV using Geometric
%%Methods on SE(3)-> no integral error terms
%Feb 19,2022

function [f, M, error, calculated] ...
    = position_control(X, desired, k, param)
% [f, M, ei_dot, eI_dot, error, calculated] = position_control(X, desired, 
% k, param)
%
% Position controller that uses decoupled-yaw controller as the attitude
% controller
% 
%   Caluclates the force and moments required for a UAV to reach a given 
%   set of desired position commands using a decoupled-yaw controller
%   defined in https://ieeexplore.ieee.org/document/8815189.
%   
%   Inputs:
%    X: (24x1 matrix) states of the system (x, v, R, W, ei, eI)
%    desired: (struct) desired states
%    k: (struct) control gains
%    param: (struct) parameters such as m, g, J in a struct
%
%  Outputs:
%    f: (scalar) required motor force
%    M: (3x1 matrix) control moment required to reach desired conditions
%    error: (struct) errors for attitude and position control (for data
%    logging)
%    calculated: (struct) calculated desired commands (for data logging)


% Unpack states
[x, v, R, W] = split_to_states(X);

m = param.m; %mass
g = param.g; %gravity
e3 = [0, 0, 1]'; %z-axis of world frame
e2 = [0, 1, 0]'; %y-axis of world frame
e1 = [1, 0, 0]'; %x-axis of world frame


error.x = x - desired.x; %eqn(17), Geometric Methods on SE(3)
error.v = v - desired.v; %eqn(18), Geometric Methods on SE(3)

A = -k.x * error.x - k.v * error.v -m * g * e3 +m * desired.x_2dot; %Appendix F, Control of Complex Maneuvers...
b3 = R * e3; % eqn(3), Decoupled Yaw Control
f = -dot(A, b3); %eqn(32), Decoupled Yaw Control
ea = g * e3 ... % eqn(41), Decoupled Yaw Control
    - f / m * b3 ...
    - desired.x_2dot ...
    + param.x_delta / m;

A_dot = - k.x * error.v ...  % derivative of A
    - k.v * ea ...
    + m * desired.x_3dot ...;

b3_dot = R * hat(W) * e3; %eqn 22 (Decoupled Yaw Control)
f_dot = -dot(A_dot, b3) - dot(A, b3_dot);
eb = - f_dot / m * b3 - f / m * b3_dot - desired.x_3dot; %time derivative of ea
A_ddot = - k.x * ea ...
    - k.v * eb ...
    + m * desired.x_4dot ...;

[b3c,b3c_dot,b3c_ddot] = deriv_unit_vector(-A,-A_dot,-A_ddot); %Appendix F, Control of Complex Maneuvers V3...

C = hat(b3c) * desired.b1; %Appendix F, Control of Complex Maaneuvers...
C_dot = hat(b3c_dot) * desired.b1 + hat(b3c) * desired.b1_dot; 
C_ddot = hat(b3c_ddot) * desired.b1 + hat(b3c_dot) * desired.b1_dot + hat(b3c_dot) * desired.b1_dot + hat(b3c) * desired.b1_2dot;

[b2c,b2c_dot,b2c_ddot] = deriv_unit_vector(-C,-C_dot,-C_ddot); % %Appendix F, Control of Complex Maneuvers...

b1c = hat(b2c) * b3c; %eqn(36), Control of Complex Maneuvers
b1c_dot = hat(b2c_dot) * b3c + hat(b2c) * b3c_dot; %Appendix F, Control of Complex Maneuvers
b1c_ddot = hat(b2c_ddot) * b3c ...
    + 2 * hat(b2c_dot) * b3c_dot ...
    + hat(b2c) * b3c_ddot;


Rc = [b1c, b2c, b3c]; %eqn (22), ...Complex Maneuvers
Rc_dot = [b1c_dot, b2c_dot, b3c_dot];
Rc_ddot = [b1c_ddot, b2c_ddot, b3c_ddot];

Wc = vee(Rc' * Rc_dot); %eqn(22), ...Complex Maneuvers... -> vee map is the inverse of hat map
Wc_dot = vee( Rc' * Rc_ddot - hat(Wc)^2); % eqn(97) ... Complex Maneuvers

eR = 1 / 2 * vee(Rc' * R - R' * Rc); %eqn(21)  ...Complex Maneuvers
eW = W - R' * Rc * Wc; %eqn(21)  ...Complex Maneuvers

M = - k.R * eR ...
    - k.W * eW ...
    + hat(R' * Rc * Wc) * param.J * R' * Rc * Wc ...
    + param.J * R' * Rc * Wc_dot; %eqn(20)  ...Complex Maneuvers


%% run attitude controller 

[M,  error.R, error.W] = attitude_control( ...
        R, W,  ...
        Rc, Wc, Wc_dot, ...
        k, param);
error.y = 0;
error.Wy = 0;

%% Saving data
calculated.b3 = b3c;
calculated.b3_dot = b3c_dot;
calculated.b3_ddot = b3c_ddot;
calculated.b1 = b1c;
calculated.R = Rc;
calculated.W = Wc;
calculated.W_dot = Wc_dot;

end