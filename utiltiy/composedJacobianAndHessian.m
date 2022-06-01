function [f, J, H] = composedJacobianAndHessian( f1, f2, x)
    % Computes the jacobian J and hess H of the function composition
    % f1(f2(x)) at the point x, by using the jacobian and hessian of f1 and
    % f2
    
    if nargout == 1
        [y] = f2(x);
        [f] = f1(y);     
    elseif nargout == 2
        [y, J2] = f2(x);
        [f, J1] = f1(y);
        J = J1*J2;
    else
        [y, J2, H2] = f2(x);
        [f, J1, H1] = f1(y);

        n = length(x);

        % Compute the pesky tensor contraction terms
        H_pesky = J1*reshape( H2, [size(H2,1) n*n]);
        H_pesky = reshape(H_pesky, [n n]); 
        
        J = J1*J2;
        H = J2'*H1*J2 + H_pesky;
    end
end