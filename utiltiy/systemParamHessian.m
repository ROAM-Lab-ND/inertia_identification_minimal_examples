function H = systemParamHessian(x, x0)
    % Hessian tensor for the systemParamMap function
    % The indicies in this function are specific to the system used
    
    [pi, b , bc] = splitParams(x);
    pi = reshape(pi,10,6);
    if nargin == 2
        pi0 = reshape( x0(1:60), 10,6);
    end
    H = zeros(66,66,66);
    for i = 1:6
        i1 = (i-1)*10+1;
        i2 = i*10;
        inds = i1:i2;
        if nargin == 1
            prior = pinertiaToVec(eye(4));
        else
            prior = pi0(:,i);
        end
        H(inds, inds, inds) = singleBodyParamHessian( pi(:,i), prior);
    end
    for i = 61:66
        H(i,i,i) = exp( x(i) );
    end
end

    