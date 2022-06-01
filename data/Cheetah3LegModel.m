function [model, graphics] = Cheetah3LegModel()

% This function creates both the rigid body model struct and a graphics
% cell array that specifies the cheetah 3 model

% The 0-configuration for the robot is with legs stright down, cheetah
% pointed along the +x axis of the ICS. 
% The inertial coordinates have +z up (i.e. gravity = [0 0 -9.81 m/s^2])
% The body coordinates have +x forward, +y left, +z up
%

model.gc.point = zeros(3,0);
model.gc.body = zeros(1,0);

Nb = 6;

%% Nominal Paramters of Cheetah
bodyWidth = .256;
bodyHeight = .200;
bodyLength = .600;


motorRad   = .06;
legRad     = 0.03;

lo = 0.045;
lc = bodyWidth/2;
l1 = .342; % Length of top link
l2 = .345; % Length of bottom link
motorWidth = 2*lo;
                  

%% Legs
side_sign = 1;

Nb = 0;
side_sign = 1;
%% Legs

%% Ab/Ad (Hip Roll)
Nb = Nb+1;

model.parent(Nb) = 0;
model.jtype{Nb}  = 'Rx';
model.Xtree{Nb}  = plux(eye(3), [0 0 0]);
model.Xrotor{Nb} = plux(eye(3), [0 -side_sign*lc 0]');

model.I{Nb}      = mcI(1.5, [0 side_sign*lo 0], boxInertia(1.5, ones(3,1)*motorWidth));
model.Irot{Nb}   = mcI(.5, [0 0 0], boxInertia(.5, [.25 1 1]'*motorWidth));

graphics{Nb}.boundCenter = [0 side_sign*lo 0]';
graphics{Nb}.boundAxes   = [motorWidth motorWidth motorWidth]/2*1.8;

graphics{Nb}.boundCenterRot = [0 0 0]';
graphics{Nb}.boundAxesRot   = [motorWidth motorWidth motorWidth]/2*1.8;

model.gr{Nb}     = 10.604;


%% Hip Pitch
Nb = Nb+1;
model.parent(Nb) = Nb-1;
model.jtype{Nb}  = 'Ry';
model.Xtree{Nb}  = plux(rz(pi), [0 side_sign*lo 0]');
model.Xrotor{Nb} = plux(rz(pi), [0 side_sign*lo 0]');

model.I{Nb}      = mcI(.5, [0 0 -l1]/2, boxInertia(.5,[legRad*2 legRad*2 l1]) );
model.Irot{Nb}   = mcI(.5, [0 0 0]/2, boxInertia(.5, [1 .25 1]'*motorWidth) );

model.gr{Nb} = 10.604;

graphics{Nb}.boundCenter = [0 0 -l1]'/2;
graphics{Nb}.boundAxes   = [legRad*2.2 legRad*2.2 l1*1.2]/2*1.8;

graphics{Nb}.boundCenterRot = [0 0 0]'/2;
graphics{Nb}.boundAxesRot   = [motorWidth motorWidth motorWidth]/2*1.8;
    
    

%% Knee Pitch
Nb = Nb+1;
model.parent(Nb) = Nb-1;
model.jtype{Nb}  = 'Ry';

model.Xtree{Nb}  = plux(eye(3), [0 0 -l1]');
model.Xrotor{Nb}  = plux(eye(3), [0 0  0 ]');

model.I{Nb}     = mcI(.5, [0 0 -l2/2], boxInertia(.5,[legRad*2 legRad*2 l2]) );
model.Irot{Nb}  = mcI(.5, [0 0 0    ], boxInertia(.5,[1 .25 1]'*motorWidth)  );
model.gr{Nb} = 10.604; 

model.NB   = Nb;

graphics{Nb}.boundCenter = [0 0 -l2]'/2;
graphics{Nb}.boundAxes   = [legRad*2 legRad*2 l2*1.2]/2*1.4;


graphics{Nb}.boundCenterRot = [0 0 0]';
graphics{Nb}.boundAxesRot   = [motorWidth motorWidth motorWidth]/2*1.8;





mT = 0;
for i = 1:length(model.I)
    mT = mT+model.I{i}(6,6);
end
%mT

model.NB   = Nb;

end

function I = boxInertia(mass, x)
    I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end