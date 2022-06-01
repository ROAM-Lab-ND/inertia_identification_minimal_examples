function [pi, Jac, Hess]= logCholeskyVecToParams(x, pi0)
    % Converts log-cholesky params to standard inertial params
    % see: Smooth Parameterization of Rigid-Body Inertia by (Rucker and
    % Wensing) https://ieeexplore.ieee.org/document/9690029
    
    d1 = x(1);
    d2 = x(2);
    d3 = x(3);
    s12 = x(4);
    s13 = x(5);
    s23 = x(6);
    t1 = x(7);
    t2 = x(8);
    t3 = x(9);
    a = x(10);
    
    U = exp(a)*[exp(d1) s12 s13, t1;0 exp(d2) s23, t2;0 0 exp(d3), t3; 0,0,0,1];
    if nargin == 1
        J0 = eye(4);
    else
        J0 = inertiaVecToPinertia(pi0);
    end
    
    J = U*J0*U.';
    
    pi = pinertiaToVec(J);
    
    if nargout > 1
        if nargin == 1
            Jac = singleBodyParamJacobian(x);
        else
            Jac = singleBodyParamJacobian(x, pi0);
        end
    end
    if nargout > 2
        if nargin == 1
            Hess = singleBodyParamHessian(x);
        else
            Hess = singleBodyParamHessian(x, pi0);
        end
    end
end

    