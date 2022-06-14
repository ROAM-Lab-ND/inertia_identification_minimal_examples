function plot_bounding_ellipsoid(T_0i, Q, trans, color)

    [Axes, center] = QtoSemiAxes(Q);
    Q = Axes*Axes';
    
    [R,S] = eig(Q);
    if det(R)<0
        R_temp = R; S_temp = S;
        R(:,1) = R_temp(:,2);
        R(:,2) = R_temp(:,1);
        S(1,1) = S_temp(2,2);
        S(2,2) = S_temp(1,1);
    end
    
    T = [R, zeros(3,1);0,0,0,1];
    T(1:3,4) = center;
    
    T = T_0i*T;
    
    S = sqrt(S);
    
    a = S(1,1); 
    b = S(2,2);
    c = S(3,3);
    
    scale = 1;
    [xc,yc,zc] = ellipsoid(0,0,0,a*scale,b*scale,c*scale,200);

    % rotate data with orientation matrix U and center M
    RR = T(1:3,1:3); M = T(1:3,4);
    
    a = kron(RR(:,1),xc); b = kron(RR(:,2),yc); c = kron(RR(:,3),zc);
    data = a+b+c; n = size(data,2);
    x = data(1:n,:)+M(1); y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
    
    sc = surf(x,y,z,'EdgeColor','none','FaceColor',color,'FaceAlpha',trans);
    