function [ G ] = calculateGecon( C, R, L, Rtrunc, Ltrunc )

    %%%%%%%%%% Now try to reproduce G in more economical manner. This
    %%%%%%%%%% doesn't work!
    
    [~, D] = size(C);

    E = C.*conj(C);
    
    lambda = diag(L(1,:))*(1./E)*diag(R(1,:));
    
    lambda_trunc = lambda(2:end, :);
    
    [Y, Dmat, Z ] = svd(lambda_trunc);
    
    Dsquare = Dmat(:, 1:D-1);
    
    Dsquare_inv = inv(Dsquare);
    
    Dinv = zeros(D, D-1);
    
    Dinv(1:D-1, 1:D-1) = Dsquare_inv;
    
   
    T = lambda*inv(Z')*Dinv*inv(Y);
    
    Rtrunc_square = Rtrunc(:, 2:end);
    
    X = inv(Rtrunc_square)*Rtrunc;
    
    XT = transpose(X).*T;
    
    W = Rtrunc*conj(XT);
    
    
    %%%%%%%%%now deal construcnt G from W and other stuff...
    
    factor_l = kron(W', L');
    
    lambda_trunc_vec = reshape(transpose(lambda_trunc), D*(D-1), 1);
    
    factor_l = diag(lambda_trunc_vec)*factor_l;
    
    factor_r = kron(transpose(Rtrunc_square), transpose(L));
    
    G = transpose(factor_l)*factor_r;
    
    
    
end

