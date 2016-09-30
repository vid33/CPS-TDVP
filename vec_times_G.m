function [ out ] = vec_times_G( in, C, L, R, Rtrunc )

    [~, D] = size(C);

    EE = C.*conj(C);

    in = transpose(reshape(in', D, D-1));
    

    T1 = in*bsxfun(@times, conj(L), R(1,:));
    T2 = transpose(conj(Rtrunc))*T1;
    T3 = T2.*(1./EE);

    out = transpose(bsxfun(@times, transpose(L(1,:)), transpose(Rtrunc)))*T3*transpose(L);

    %out = 
    
    out = reshape(transpose(out), 1, D*(D-1));

end

