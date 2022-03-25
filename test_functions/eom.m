function Xdot = eom(t, X, k, param)

e3 = [0, 0, 1]';
m = param.m;
J = param.J;

[~, v, R, W] = split_to_states(X);

A = 1;
B = 1;
C = 0.2;

d = pi / 2 * 0;

a = 1;
b = 2;
c = 2;
alt = -1;
m = 0:0.05:10;

waypts = [A * sin(a *m + d);B * sin(b * m);alt + C * cos(2 * m)];

desired = command_lissajous(t,waypts);
    
[f, M, ~, ~] = position_control(X, desired, k, param);

xdot = v;
vdot = param.g * e3 ...
    - f / m * R * e3 + param.x_delta / m;
Wdot = J \ (-hat(W) * J * W + M + param.R_delta);
Rdot = R * hat(W);

Xdot=[xdot; vdot; Wdot; reshape(Rdot,9,1)];
end