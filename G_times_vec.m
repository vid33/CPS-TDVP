function [ out ] = G_times_vec( in, C, L, R, Rtrunc )

    [~, D] = size(C);

    EE = C.*conj(C);
    
    in = transpose(reshape(in, D, D-1));

    %in = transpose(reshape(in, D, D-1));
    
    T1 = transpose(in)*transpose(bsxfun(@times, transpose(L(1,:)), transpose(Rtrunc)));
    T2 = transpose(T1)*L;
    T3 = T2.*(1./EE);

    out = conj(Rtrunc)*T3*transpose(bsxfun(@times, conj(L), R(1,:)));
    
    out = reshape(transpose(out), D*(D-1), 1);

end

