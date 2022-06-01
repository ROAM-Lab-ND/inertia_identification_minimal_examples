function [f, J, H] = regularizedLeastSquares( x, data)
    % Computes the objective for an entropically regularized least-square
    % problem
    % The input x represents the standard params. 
    % Note: The indicies in this function are specific to the system used
    
    Y = data.Y;
    e = Y*x - data.tau ;
    
    if nargout > 1
        J = e.'*Y / length(e) ;
    end
    if nargout > 2
        H = data.Y'*data.Y / length(e) ;
    end
    
    f = 1/2*e.'*e./length(e);
    pi = reshape(x(1:60), 10,6); % Specific to this system
    
    for i = 1:6
       i1 = (i-1)*10+1;
       i2 = i*10;
       inds = i1:i2;
       if nargout == 1
          f = f+ data.gamma*entropicDivergence( pi(:,i), data.prior(:,i) );
       elseif nargout == 2
          [f_reg, J_reg] = entropicDivergence( pi(:,i), data.prior(:,i) );
          
          f = f + data.gamma*f_reg;
          J(inds) = J(inds) + data.gamma*J_reg.';
       else
          [f_reg, J_reg, H_reg] = entropicDivergence( pi(:,i), data.prior(:,i) );
          f = f + data.gamma*f_reg;
          J(inds) = J(inds) +  data.gamma*J_reg.';
          H(inds, inds) = H(inds, inds) + data.gamma*H_reg;
       end    
    end
end
    