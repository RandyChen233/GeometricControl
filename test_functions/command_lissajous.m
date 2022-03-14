
function desired = command_lissajous(t)
addpath('C:\Users\Randy666\Documents\projects\AdvancedControlsLab\L1 adaptive+geometric\minimum_snap');
%% condition
A = 1;
B = 1;
C = 0.2;

d = pi / 2 * 0;

a = 1;
b = 2;
c = 2;
alt = -1;
t = 0:0.1:25; %time
waypts = [A * sin(a *t + d);B * sin(b * t);alt + C * cos(2 * t)];

v0 = [0,0,0];
a0 = [0,0,0];
v1 = [0,0,0];
a1 = [0,0,0];
T = t(end);
ts = arrangeT(waypts,T);
n_order = 7;
% trajectory plan
polys_x = minimum_snap_single_axis_simple(waypts(1,:),ts,n_order,v0(1),a0(1),v1(1),a1(1));
polys_y = minimum_snap_single_axis_simple(waypts(2,:),ts,n_order,v0(2),a0(2),v1(2),a1(2));
polys_z = minimum_snap_single_axis_simple(waypts(3,:),ts,n_order,v0(3),a0(3),v1(3),a1(3));

xx = [];
yy = [];
zz = [];
vxx = [];
vyy = [];
vzz = [];
axx = [];
ayy = [];
azz = [];
jxx = [];
jyy = [];
jzz = [];
sxx = [];
syy = [];
szz = [];


for i=1:size(polys_x,2)
    tt = ts(i):0.01:ts(i+1);
    xx = [xx,polys_vals(polys_x,ts,tt,0)];
    yy = [yy,polys_vals(polys_y,ts,tt,0)];
    zz = [zz,polys_vals(polys_z,ts,tt,0)];

    vxx = [vxx,polys_vals(polys_x,ts,tt,1)];
    axx = [axx,polys_vals(polys_x,ts,tt,2)];
    jxx = [jxx,polys_vals(polys_x,ts,tt,3)];
    sxx = [sxx,polys_vals(polys_x,ts,tt,4)];

    vyy = [vyy,polys_vals(polys_y,ts,tt,1)];
    ayy = [ayy,polys_vals(polys_y,ts,tt,2)];
    jyy = [jyy,polys_vals(polys_y,ts,tt,3)];
    syy = [syy,polys_vals(polys_y,ts,tt,4)];

    vzz = [vzz,polys_vals(polys_z,ts,tt,1)];
    azz = [azz,polys_vals(polys_z,ts,tt,2)];
    jzz = [jzz,polys_vals(polys_z,ts,tt,3)];
    szz = [szz,polys_vals(polys_z,ts,tt,4)];
    
end



%% saving data

desired.x = [xx,yy,zz]';
desired.v = [vxx,vyy,vzz]';
desired.x_2dot = [axx,ayy,azz]';
desired.x_3dot = [jxx,jyy,jzz]';
desired.x_4dot = [sxx,syy,szz]';
% 
% desired.x = [A * sin(a *t + d),B * sin(b * t),alt + C * cos(2 * t)]';
% 
% desired.v = [A * a * cos(a * t + d), ...
%    B * b * cos(b * t), ...
%    C * c * -sin(c * t)]';
% 
% desired.x_2dot = [A * a^2 * -sin(a * t + d), ...
%    B * b^2 * -sin(b * t), ...
%    C * c^2 * -cos(c * t)]';
% 
% desired.x_3dot = [A * a^3 * -cos(a * t + d), ...
%    B * b^3 * -cos(b * t), ...
%    C * c^3 * sin(c * t)]';
% 
% desired.x_4dot = [A * a^4 * sin(a * t + d), ...
%     B * b^4 * sin(b * t), ...
%     C * c^4 * cos(c * t)]';

w = 2 * pi / 10;
%t = linspace(0,t(end),size(xx,2));
desired.b1 = [cos(w * t),sin(w * t), zeros(1,length(t))]';
desired.b1_dot = w * [-sin(w * t),cos(w * t),zeros(1,length(t))]';
desired.b1_2dot = w^2 * [-cos(w * t),-sin(w * t), zeros(1,length(t))]';


end