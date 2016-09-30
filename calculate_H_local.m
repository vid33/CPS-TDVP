function [ Hlocal  ] = calculate_H_local( siteOverlap, sx, sz, J, h )

    Hlocal = cell(4, 1);
    
    Hlocal{1} = zeros(2^siteOverlap);
    
    for kk=1:siteOverlap
        Hlocal{1} =  Hlocal{1} + (h/(2*siteOverlap))*kron(kron(eye(2^(kk-1)), sz), eye(2^(siteOverlap-kk))) ;
    end
    
    for kk=1:(siteOverlap-1)
        Hlocal{1} = Hlocal{1} - (J/(2*siteOverlap))*kron(kron(eye(2^(kk-1)), kron(sx,sx)), eye(2^(siteOverlap-1-kk)));
    end
   
    Hlocal{4} = Hlocal{1};
    
    %split ones
    Hlocal{2} = (sqrt(J/(siteOverlap)))*kron(eye(2^(siteOverlap-1)), sx);
    Hlocal{3} =  -(sqrt(J/(siteOverlap)))*kron(sx, eye(2^(siteOverlap-1)));

end

