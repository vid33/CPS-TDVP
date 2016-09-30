function [ twoPointLocal  ] = calculate_twoPoint_local( siteOverlap, sz, PP)
    %two-point sz  - - PP - - - sz operator  PP>3 PP<siteOverlap)

    twoPointLocal = cell(1 + (PP-1)*2, 1);
    
    twoPointLocal{end} = zeros(2^siteOverlap);
    
    for kk=1:(siteOverlap - (PP-1) )
        tmp = (1/(2*siteOverlap))*kron(kron(eye(2^(kk-1)), sz), eye(2^(PP-2)));
        tmp = kron(tmp, sz);
        twoPointLocal{end} = twoPointLocal{end} + kron(tmp, eye(2^(siteOverlap-(kk-1) - PP)) );
        %twoPointLocal{end} =  twoPointLocal{end} + kron(kron(eye(2^(kk-1)), sz), eye(2^(siteOverlap-kk))) ;
    end
    
   % for kk=1:(siteOverlap-1)
   %     Hlocal{1} = Hlocal{1} - (J/(2*siteOverlap))*kron(kron(eye(2^(kk-1)), kron(sx,sx)), eye(2^(siteOverlap-1-kk)));
   % end
   
    for kk=1:(PP-1)
        twoPointLocal{kk} = (sqrt(1/(siteOverlap)))*kron(eye(2^(siteOverlap-1-(kk-1))), sz);
        twoPointLocal{kk} = kron(twoPointLocal{kk}, eye(2^(kk-1)));
        
        twoPointLocal{PP-1+kk} = (sqrt(1/(siteOverlap)))*kron(eye(2^(kk-1)), sz);
        twoPointLocal{PP-1+kk} = kron(twoPointLocal{PP-1+kk}, eye(2^(siteOverlap-1-(kk-1))));   
    end
    
    %Hlocal{4} = Hlocal{1};
    
    %split ones
    %Hlocal{2} = (sqrt(1/(siteOverlap)))*kron(eye(2^(siteOverlap-1)), sx);
    %Hlocal{3} =  -(sqrt(1/(siteOverlap)))*kron(sx, eye(2^(siteOverlap-1)));

end

