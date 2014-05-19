function P = kPermZmZn(M,N)

P = zeros(2,M*N);
A = randperm(M*N);

for j = 1:M*N

    x = A(j);

    P(1,j) = kMod(x,M);
    P(2,j) = ceil(x/M);

end