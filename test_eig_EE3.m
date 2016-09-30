bigDelta = zeros(D, D, D);

A = cpxrand(D,D);
B = cpxrand(D,D);


for kk=1:D
    bigDelta(kk,kk,kk)=1;
end

bigC = Contract({C, bigDelta}, {[-1, 1], [1, -2, -3]});

bigCbar = Contract({conj(C), bigDelta}, {[-1, 1], [1 -2, -3]});

bigEE = Contract({bigC, bigCbar}, {[-1, 1, -3], [-2, 1, -4]});


bigEEmat = transpose(reshape(bigEE, D^2, D^2));



bigA = Contract({A, bigDelta}, {[-1, 1], [1, -2, -3]});

bigB = Contract({B, bigDelta}, {[-1, 1], [1, -2, -3]});


bigEEmix = Contract({bigA, bigB}, {[-1, 1, -3], [-2, 1, -4]});


bigEEmix_mat = transpose(reshape(bigEEmix, D^2, D^2));


EEmix = A.*B;


