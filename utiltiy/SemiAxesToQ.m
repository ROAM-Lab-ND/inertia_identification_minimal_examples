function Q = SemiAxesToQ(Axes,center)
    if nargin == 1
        center = [0;0;0];
    end
    
    Q = inv(Axes*Axes');
    
    Qc  = Q*center;
    Q = [ -Q Qc ; Qc' 1-center'*Qc];
end

