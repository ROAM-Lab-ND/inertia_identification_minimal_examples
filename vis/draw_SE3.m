function draw_SE3(T,string_color)
p=T(1:3,4);
ax=p+T(1:3,1)*0.12;
ay=p+T(1:3,2)*0.12;
az=p+T(1:3,3)*0.12;
hold on;
if nargin < 2
plot3([p(1),ax(1)], [p(2),ax(2)], [p(3),ax(3)],'black','LineWidth', 1.5);
plot3([p(1),ay(1)], [p(2),ay(2)], [p(3),ay(3)],'black','LineWidth', 1.5);
plot3([p(1),az(1)], [p(2),az(2)], [p(3),az(3)],'black','LineWidth', 1.5);
else
plot3([p(1),ax(1)], [p(2),ax(2)], [p(3),ax(3)],string_color,'LineWidth', 1.5);
plot3([p(1),ay(1)], [p(2),ay(2)], [p(3),ay(3)],string_color,'LineWidth', 1.5);
plot3([p(1),az(1)], [p(2),az(2)], [p(3),az(3)],string_color,'LineWidth', 1.5);    
end
% axis([-1/2 1/2 -1/2 1.5/2 0 1.2]);
end