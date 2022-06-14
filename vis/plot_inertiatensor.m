function sc = plot_inertiatensor(T_sb_i, G_b_i, trans, color, varargin)
n_parts = size(G_b_i,3);

for i=[1:n_parts]
    % [R,S,V] = svd(G_sbar_init{i}(1:3,1:3));
    % a = sqrt(5*(-S(1,1)+S(2,2)+S(3,3))/(2*G_sbar_init{i}(4,4)));
    % b = sqrt(5*(S(1,1)-S(2,2)+S(3,3))/(2*G_sbar_init{i}(4,4)));
    % c = sqrt(5*(S(1,1)+S(2,2)-S(3,3))/(2*G_sbar_init{i}(4,4)));
    T_bc = [eye(3,3), skew(G_b_i(1:3,4:6,i))./G_b_i(4,4,i);0,0,0,1];
    G_c = AdjointMatrix(T_bc)'*G_b_i(:,:,i)*AdjointMatrix(T_bc);
    % [R,S,V] = svd(G_b(1:3,1:3));
    [R,S,V] = eig(G_c(1:3,1:3));
    if det(R)<0
        R_temp = R; S_temp = S;
        R(:,1) = R_temp(:,2);
        R(:,2) = R_temp(:,1);
        S(1,1) = S_temp(2,2);
        S(2,2) = S_temp(1,1);
    end
    % 0.5*eye(3,3)*trace(S) - S
    T = [R, zeros(3,1);0,0,0,1];
    T(1:3,4) = skew(G_b_i(1:3,4:6,i))./G_b_i(4,4,i);
    T = T_sb_i(:,:,i)*T;
    a = sqrt(5*(-S(1,1)+S(2,2)+S(3,3))/(2*(G_b_i(4,4,i))));
    if ~isreal(a)
        a = 0
    end
    b = sqrt(5*(S(1,1)-S(2,2)+S(3,3))/(2*(G_b_i(4,4,i))));
    if ~isreal(b)
        b = 0
    end
    c = sqrt(5*(S(1,1)+S(2,2)-S(3,3))/(2*(G_b_i(4,4,i))));
    if ~isreal(c)
        c = 0
    end
    % T(1:3,1:3) = T(1:3,1:3)*[a,0,0;0,b,0;0,0,c];
    % draw_SE3(T);
    % transl = T(1:3,1:3)'*T(1:3,4);
    % density = G_b_i(4,4,i)/((a+0.001)*(b+0.001)*(c+0.001)*2e3);
    scale = 1;
    [xc,yc,zc] = ellipsoid(0,0,0,a*scale,b*scale,c*scale,200);

    % rotate data with orientation matrix U and center M
    RR = T(1:3,1:3); M = T(1:3,4);
    a = kron(RR(:,1),xc); b = kron(RR(:,2),yc); c = kron(RR(:,3),zc);
    data = a+b+c; n = size(data,2);
    x = data(1:n,:)+M(1); y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
    % scatter3(M(1),M(2),M(3), 'MarkerEdgeColor',color(i,:),'MarkerFaceColor',color(i,:));
    % now plot the rotated ellipse
    % color_z = ones(size(z))*color;

    % text(M(1),M(2),M(3),[ num2str(G_b_i(4,4,i)) 'kg'],'HorizontalAlignment','left','FontSize',8);

    if nargin > 4
        sc_in = varargin{1};
        set(sc_in,'XData' ,x);
        set(sc_in,'YData' ,y);
        set(sc_in,'ZData' ,z);
    else
        sc = surf(x,y,z,'EdgeColor','none','FaceColor',color(i,:),'FaceAlpha',trans);
    end


    end
    % shading interp;
    % lighting phong;
    % camlight;
    axis equal; 
    % drawnow; %hold off;
    % draw_SE3(T_sb_i);
end

function Ad_T = AdjointMatrix(T)
    R = T(1:3,1:3);
    p = T(1:3,4);
    Ad_T = [R, zeros(3,3);skew(p)*R, R];
end