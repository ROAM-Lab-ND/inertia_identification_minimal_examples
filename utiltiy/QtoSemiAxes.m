function [Axes, center] = QtoSemiAxes(Q)
    [Axes, D] = eig(-Q(1:3,1:3));
    Axes = Axes/sqrt(D);
    Qc = Q(1:3,4);
    Q  = inv(-Q(1:3,1:3));
    center = Q*Qc;
end

