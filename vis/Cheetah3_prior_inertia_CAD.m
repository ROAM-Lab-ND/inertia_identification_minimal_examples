function [pi_prior, P, Q_bound_ellipse] = Cheetah3_prior_inertia_CAD()

    [model, graphics]     = Cheetah3LegModel();
    NB = model.NB;
    
    pi_prior = zeros(10, NB*2);
    Q_bound_ellipse = cell(NB*2,1);
    P =cell(NB*2,1);    % Pseudo inertia matrix for body i
    Qc=cell(NB*2,1);    % 3x3 Matrix which provides a bounding ellipse for the full rigid body
    c =cell(NB*2,1);    % Center of Bounding Ellipsoid

    % Prior Params
    for i = 1 : model.NB
       pi_prior(:,i)   = inertiaMatToVec(model.I{i});
       pi_prior(:,i+3) = inertiaMatToVec(model.I_rotor{i});
    end

    % Loop through one body at a time.
    for i = 1:2*NB
        P{i} = inertiaVecToPinertia(pi_prior(:,i));

        assert( min(eig(P{i})) > 0 )
        
        if i <= NB
            boundSemiAxisLengths   = graphics{i}.boundAxes;
            boundCenter            = graphics{i}.boundCenter;
        else
            boundSemiAxisLengths   = graphics{i-NB}.boundAxesRot;
            boundCenter            = graphics{i-NB}.boundCenterRot;
        end

        Q_bound_ellipse{i} = SemiAxesToQ( diag(boundSemiAxisLengths) ,boundCenter);

        % If the parameters a{i} are realizable with density on the bounding
        % ellipse, this trace must be positive
        tQP = trace(Q_bound_ellipse{i}*P{i});
        assert( tQP > 0 )
    end
end
