function [ sxlocal  ] = calculate_sx_local( siteOverlap, sx )

    sxlocal = cell(4, 1);
    
    sxlocal{1} = zeros(2^siteOverlap);
    
    for kk=1:siteOverlap
        sxlocal{1} =  sxlocal{1} + (1/(2*siteOverlap))*kron(kron(eye(2^(kk-1)), sx), eye(2^(siteOverlap-kk))) ;
    end
    
   % for kk=1:(siteOverlap-1)
    %    Hlocal{1} = Hlocal{1} - (J/(2*siteOverlap))*kron(kron(eye(2^(kk-1)), kron(sx,sx)), eye(2^(siteOverlap-1-kk)));
    %end
   
    sxlocal{4} = sxlocal{1};
    
    %split ones
    sxlocal{2} = (sqrt(1/(siteOverlap)))*kron(eye(2^(siteOverlap-1)), zeros(2));
    sxlocal{3} =  -(sqrt(1/(siteOverlap)))*kron(zeros(2), eye(2^(siteOverlap-1)));

end

