function [y, J, H] = systemParamMap(x, x0)
    % This function maps from log cholesky params to conventional params
    % The indicies in this function are specific to the system used

    [pi, b , bc] = splitParams(x);
    pi = reshape(pi,10,6);

    if nargin == 2
        pi0 = reshape( x0(1:60), 10,6);
    end
    
    for i = 1:6
        if nargin == 1
            y(:,i) = logCholeskyVecToParams( pi(:,i) );
        else
            y(:,i) = logCholeskyVecToParams( pi(:,i), pi0(:,i) );
        end  
    end
    y = y(:);
    
    % Friction coefficients
    for i = 61:66
        y(i) = exp( x(i) );
    end
    
    if nargout > 1
        J = systemParamJacobian(x, x0);
    end
    if nargout >2
        H = systemParamHessian(x, x0);
    end
    
end