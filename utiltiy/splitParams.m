function [pi, b , bc] = splitParams(x)
   % Breaks out parameters into those of rigid bodies and for
   % viscous/coulomb friction
   % Note: The indicies in this function are specific to the system used
 
   pi = x(1:60); % Rigid body params
   b  = x(61:63); % viscous friction
   bc = x(64:66); % coulomb friction
end
