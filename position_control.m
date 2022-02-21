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

error.x = x - desired.x; %eqn(17)
error.v = v - desired.v; %eqn(18)


f = dot((k.x * error.x + k.v * error.v + m * g * e3 - m * desired.x_2dot),R * e3); %eqn(19)

A = dot(-f,R*e3); % Section III. and eqn(14) in Geometric Tracking Control.....

%I just realized that Control of Complex Maneuvers for a Quadrotor UAV using Geometric
%Methods on SE(3) ,  relies one some results obtained in Geometric Tracking
%Control(2010).... I have to reference some results in Geometric Tracking Control
%to obtain the time derivatives of Rc, the computed attitude ; otherwise, W_c
%cannot be established and hence cannot be fed into attitude_control()
%However...if I use this ficticious control input representation A
%according to Geometric Controls of a Quadrotor UAV with Decoupled Yaw
%Control, all the changes I've made is useless because I only needed to
%remove the integral action in the control signal A at the very
%beginning??????what

[b3c] = deriv_unit_vector(-(k.x * error.x + k.v * error.v + m * g * e3 - m * desired.x_2dot)); %eqn(23)
b1c = -cross(b3c,cross(b3c,desired.x))/norm(cross(b3c,desired.x)); %eqn(38)
b2c = cross(b3c,b1c); 
Rc = [b1c, b2c, b3c]; %eqn(22)

M = - k.R * eR ...
    - k.W * eW ...
    + hat(R' * Rc * Wd) * param.J * R' * Rc * Wd ...
    + param.J * R' * Rc * Wd_dot;  %eqn(20)

Wc = vee(Rc' * Rc_dot);%eqn(22) -> vee map is the inverse of hat map
%Wc_dot = vee(Rc' * Rc_ddot - hat(Wc)^2); %1st derivative of Wc


%% run attitude controller 
[f,M, error.R, error.W] = attitude_control( ...
    X,R, W,  ...  % states
    Rc, Wc, ...  % desired values
    k, param ...  % gains and parameters
)

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
%calculated.W3 = dot(R * e3, Rc * Wc);
%calculated.W3_dot = dot(R * e3, Rc * Wc_dot) ...
%    + dot(R * hat(W) * e3, Rc * Wc);
end