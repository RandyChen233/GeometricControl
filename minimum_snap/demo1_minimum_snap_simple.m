function demo1_minimum_snap_simple()
    clear,clc;

    %% condition
    A = 1;
    B = 1;
    C = 0.2;
    
    d = pi / 2 * 0;
    
    a = 1;
    b = 2;
    c = 2;
    alt = -1;
    m = 0:0.1:25; %time
    waypts = [A * sin(a *m+ d);B * sin(b * m);alt + C * cos(2 * m)];

    v0 = [0,0,0];
    a0 = [0,0,0];
    v1 = [0,0,0];
    a1 = [0,0,0];
    T = m(end);
    ts = arrangeT(waypts,T);
    n_order = 7;
    
    %% trajectory plan
    polys_x = minimum_snap_single_axis_simple(waypts(1,:),ts,n_order,v0(1),a0(1),v1(1),a1(1));
    polys_y = minimum_snap_single_axis_simple(waypts(2,:),ts,n_order,v0(2),a0(2),v1(2),a1(2));
    polys_z = minimum_snap_single_axis_simple(waypts(3,:),ts,n_order,v0(3),a0(3),v1(3),a1(3));
    
    
    %% result show
    figure(1)
    %plot3(waypts(1,:),waypts(2,:),waypts(3,:),'*r');hold on;
    plot3(waypts(1,:),waypts(2,:),waypts(3,:),'b--');
    title('minimum snap trajectory');
    color = ['grc'];
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
    figure(2)
    
    ts = linspace(0,ts(end),size(xx,2));
    subplot(5,3,1),plot(ts,xx);title('x position');
    subplot(5,3,2),plot(ts,yy);title('y position');
    subplot(5,3,3),plot(ts,zz);title('z position');
    subplot(5,3,4),plot(ts,vxx);title('x velocity');
    subplot(5,3,5),plot(ts,vyy);title('y velocity');
    subplot(5,3,6),plot(ts,vzz);title('z velocity');
    subplot(5,3,7),plot(ts,axx);title('x acceleration');
    subplot(5,3,8),plot(ts,ayy);title('y acceleration');
    subplot(5,3,9),plot(ts,azz);title('z acceleration');
    subplot(5,3,10),plot(ts,jxx);title('x jerk');
    subplot(5,3,11),plot(ts,jyy);title('y jerk');
    subplot(5,3,12),plot(ts,jzz);title('z jerk');
    subplot(5,3,13),plot(ts,sxx);title('x snap');
    subplot(5,3,14),plot(ts,syy);title('y snap');
    subplot(5,3,15),plot(ts,szz);title('z snap');
end


function polys = minimum_snap_single_axis_simple(waypts,ts,n_order,v0,a0,ve,ae)
p0 = waypts(1);
pe = waypts(end);

n_poly = length(waypts)-1;
n_coef = n_order+1;

% compute Q
Q_all = [];
for i=1:n_poly
    Q_all = blkdiag(Q_all,computeQ(n_order,4,ts(i),ts(i+1)));
end
b_all = zeros(size(Q_all,1),1);

Aeq = zeros(4*n_poly+2,n_coef*n_poly);
beq = zeros(4*n_poly+2,1);

% start/terminal pva constraints  (6 equations)
Aeq(1:3,1:n_coef) = [calc_tvec(ts(1),n_order,0);
                     calc_tvec(ts(1),n_order,1);
                     calc_tvec(ts(1),n_order,2)];
Aeq(4:6,n_coef*(n_poly-1)+1:n_coef*n_poly) = ...
                    [calc_tvec(ts(end),n_order,0);
                     calc_tvec(ts(end),n_order,1);
                     calc_tvec(ts(end),n_order,2)];
beq(1:6,1) = [p0,v0,a0,pe,ve,ae]';

% mid p constraints    (n_ploy-1 equations)
neq = 6;
for i=1:n_poly-1
    neq=neq+1;
    Aeq(neq,n_coef*i+1:n_coef*(i+1)) = calc_tvec(ts(i+1),n_order,0);
    beq(neq) = waypts(i+1);
end

% continuous constraints  ((n_poly-1)*3 equations)
for i=1:n_poly-1
    tvec_p = calc_tvec(ts(i+1),n_order,0);
    tvec_v = calc_tvec(ts(i+1),n_order,1);
    tvec_a = calc_tvec(ts(i+1),n_order,2);
    neq=neq+1;
    Aeq(neq,n_coef*(i-1)+1:n_coef*(i+1))=[tvec_p,-tvec_p];
    neq=neq+1;
    Aeq(neq,n_coef*(i-1)+1:n_coef*(i+1))=[tvec_v,-tvec_v];
    neq=neq+1;
    Aeq(neq,n_coef*(i-1)+1:n_coef*(i+1))=[tvec_a,-tvec_a];
end

Aieq = [];
bieq = [];

p = quadprog(Q_all,b_all,Aieq,bieq,Aeq,beq);

polys = reshape(p,n_coef,n_poly);

end

