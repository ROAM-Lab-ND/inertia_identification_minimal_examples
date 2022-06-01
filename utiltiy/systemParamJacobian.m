function J = systemParamJacobian(x, x0)
    % Jacobian of the systemParamMap function
    % The indicies in this function are specific to the system used
    
    pi = reshape(x(1:60),10,6);
    if nargin == 2
        pi0 = reshape( x0(1:60), 10,6);
    end
    
    for i = 1:6
        i1 = (i-1)*10+1;
        i2 = i*10;
        inds = i1:i2;
        
        if nargin == 1
            prior = pinertiaToVec(eye(4));
        else
            prior = pi0(:,i);
        end    
        J( inds, inds ) = singleBodyParamJacobian( pi(:,i), prior);
    end
    J(61:66, 61:66) = diag( exp(x(61:end)) );
end

    