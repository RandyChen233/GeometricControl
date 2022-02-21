function desired = command_lissajous(t)

A = 1;
B = 1;
C = 0.2;

d = pi / 2 * 0;

a = 1;
b = 2;
c = 2;
alt = -1;

% t = linspace(0, 2*pi, 2*pi*100+1);
% x = A * sin(a * t + d);
% y = B * sin(b * t);
% z = alt + C * cos(2 * t);
% plot3(x, y, z);

desired.x = [A * sin(a * t + d), B * sin(b * t), alt + C * cos(2 * t)]';

desired.v = [A * a * cos(a * t + d), ...
    B * b * cos(b * t), ...
    C * c * -sin(c * t)]';

desired.x_2dot = [A * a^2 * -sin(a * t + d), ...
    B * b^2 * -sin(b * t), ...
    C * c^2 * -cos(c * t)]';

desired.x_3dot = [A * a^3 * -cos(a * t + d), ...
    B * b^3 * -cos(b * t), ...
    C * c^3 * sin(c * t)]';

desired.x_4dot = [A * a^4 * sin(a * t + d), ...
    B * b^4 * sin(b * t), ...
    C * c^4 * cos(c * t)]';

desired.w = [2 * pi / 10 , 0 , 0 ]'; %desired angular velocity
desired.w_dot = 0;
desired.w_ddot = 0;

desired.b1 = [cos(desired.w * t), sin(desired.w * t), 0]';
desired.b1_dot = desired.w * [-sin(desired.w * t), cos(desired.w * t), 0]';
desired.b1_2dot = desired.w^2 * [-cos(desired.w * t), -sin(desired.w * t), 0]';

end