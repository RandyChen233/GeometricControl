[f,M, ev, eW, eR] = attitude_control( ...
    X,R, W,  ...  % states
    Rc, Wc, Wc_dot ...  % desired values
    k, param ...  % gains and parameters
)

[x, v, R, W] = split_to_states(X);

e3 = [0, 0, 1]'; %z-axis of world frame
e2 = [0, 1, 0]'; %y-axis of world frame
e1 = [1, 0, 0]'; %x-axis of world frame

m = param.m;
g = param.g;


ev = v - desired.v; %eqn 37
eR = 1 / 2 * vee(Rc' * R - R' * Rc); %eqn 40
eW = W - R' * Rc * Wc; %eqn 40


f = dot((k.v * ev + m * g * e3 - m * desired.x_2dot),R*e3); %eqn 38


A = -k.v * ev -m * g * 3 +_ m * desired.x_2dot; %eqn 42;
A_dot = -k.v *( -desired.x2_dot) + m * desired.x_3dot; 
A_ddot =-k.v * (-desired.x2_dot) +  m * desired.x_4dot;

[b3c,b3c_dot,b3c_ddot] = deriv_unit_vector(-A,-A_dot,-A_ddot); %eqn(23)




M = - kR * eR ...
    - kW * eW ...
    + hat(R' * Rc * Wc) * param.J * R' * Rc * Wc ...
    + param.J * R' * Rc * Wc_dot;  %eqn 39



