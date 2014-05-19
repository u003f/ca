function [cellsAll,gProfileAll,cProfileAll,hProfileAll,phiaAll,...
    mutant1All,mutant2All,mutant3All] = kCA(X)

%If no input, use defualt parameters
if (nargin == 0)
    X(1) = 10; X(2) = 8.6e3; X(3) = 1/3; X(4) = 1e-3;
end

%Parameter values
k = X(1); hT = X(2); pc = X(3); pa = X(4);


nA = 1/18; dg = 130; dc = 5; a0 = 0.1; hN = 9.3e2;

N = 50;
genEnd = 5000;
mutEnd = 0.95;

hyperApprox = 1; %Approximate phase until hyperplasia?
glucApprox = 0; %Approximate glucose distribution?

% EDIT, 19 MAY 2014: pc < 1, glucApprox = 0 not allowed.
pc = 1; glucApprox = 0;
% /EDIT

%--------------------------------------------------------------------------
%Initial conditions.
[cells,mutant1,mutant2,mutant3,genNum] = kInit(hyperApprox,N,pa);

cellsAll = []; gProfileAll = []; cProfileAll = []; hProfileAll = [];
phiaAll = []; mutant1All = []; mutant2All = []; mutant3All =[];

finished = 0;

while (finished == 0)

    genNum = genNum + 1;

    cells2 = max(cells,0);

    [gProfile,cProfile,hProfile,phia] = kUpdateMetabolites(cells2,...
        mutant2,nA,k,dg,dc,glucApprox);
    
    cellsAll = kUpdateAll(cells,cellsAll);
    gProfileAll = kUpdateAll(gProfile,gProfileAll);
    cProfileAll = kUpdateAll(cProfile,cProfileAll);
    hProfileAll = kUpdateAll(hProfile,hProfileAll);
    phiaAll = kUpdateAll(phia,phiaAll);
    mutant1All = kUpdateAll(mutant1,mutant1All);
    mutant2All = kUpdateAll(mutant2,mutant2All);
    mutant3All = kUpdateAll(mutant3,mutant3All);

    if (genNum <= genEnd)

        [cells,mutant1,mutant2,mutant3] = kUpdateCells(cells,gProfile,...
            cProfile,hProfile,phia,mutant1,mutant2,mutant3,a0,hN,hT,pc,...
            pa,k,nA);

    end
    
    cells2 = max(cells,0);

    if (sum(sum(cells2)) == 0)
        R = 0;
        break
    end

    mutNum = sum(sum(mutant1.*mutant2.*mutant3))/sum(sum(cells2));

    if ((genNum >= genEnd)|(mutNum >= mutEnd))
        finished = 1;
        if (genNum >= genEnd)
            R = 0;
        else
            R = (genNum - 1)^(-1);
        end
    end
end

%--------------------------------------------------------------------------
function [cells,mutant1,mutant2,mutant3,genNum] = kInit(hyperApprox,N,pa)

cells = ones(1,N);

if (hyperApprox == 0)
    genNum = 0;
    mutant1 = zeros(1,N); mutant2 = zeros(1,N); mutant3 = zeros(1,N);

else
    p = 1 - (1 - pa)^N; %Probability of at least one cell mutating /gen
    genNum = 1 + kGeornd(p);

    mutant1 = kMutInit(1,N);
    mutant2 = kMutInit(kBinornd(genNum-1,p),N);
    mutant3 = kMutInit(kBinornd(genNum-1,p),N);
end
%--------------------------------------------------------------------------
function mutant = kMutInit(numMut,N);

mutant = zeros(1,N);

if (numMut ~= 0)
    for k = 1:numMut
        r = kRandint(1,1,[1,N]); mutant(r) = mod(mutant(r)+1,2);
    end
end
%--------------------------------------------------------------------------
function [gProfile,cProfile,hProfile,phia] = kUpdateMetabolites(...
    cells2,mutant2,n,k,dg,dc,glucApprox)

deltag = (cells2 + (k-1)*mutant2)/dg^2;
deltac = cells2/dc^2;

if (glucApprox == 1)
    gProfile = gApprox(deltag);
else
    gProfile = kProfile1(deltag);
end

cProfile = kProfile1(deltac);

phig = (cells2 + (k-1)*mutant2).*gProfile;
phic = cells2.*cProfile;
phia = phic + n*(phig - phic);
phih = phig - phic;
if (min(min(phih))<0)
    disp('error: phic > phig occurs within model')
end

hProfile = kProfile2(phih);
%--------------------------------------------------------------------------
function XProfile = kProfile1(deltaX)

[M,N] = size(deltaX);

A = kMakeMatrix(M,N);

for m = 1:M
    for n = 1:N
        deltaX2(N*(m-1) + n) = deltaX(m,n);
    end
end

A = sparse(A - diag(deltaX2));

b = sparse(M*N,1);
b(1:N) = -1;

X0 = A\b;

for m = 1:M
    for n = 1:N
        XProfile(m,n) = X0(N*(m-1) + n);
    end
end
%--------------------------------------------------------------------------
function A = kMakeMatrix(M,N)

A0 = -4*speye(N);
A0(1,N) = 1;
A0(N,1) = 1;
for k = 1:(N-1)
    A0(k,k+1) = 1;
    A0(k+1,k) = 1;
end

I = speye(N);

A = sparse(M*N,M*N);

for m = 1:(M-1)
    K = ((m-1)*N + 1):m*N;
    Kp = (m*N + 1):(m+1)*N;
    A(K,K) = A0;
    A(K,Kp) = I;
    A(Kp,K) = I;
end

KM = ((M-1)*N + 1):M*N;

A(KM,KM) = A0 + I;
%--------------------------------------------------------------------------
function hProfile = kProfile2(phih)

[M,N] = size(phih);

A = kMakeMatrix(M,N);

phih2 = zeros(M*N,1);

for m = 1:M
    for n = 1:N
        phih2(N*(m-1) + n) = -phih(m,n);
    end
end

h0 = A\phih2;

for m = 1:M
    for n = 1:N
        hProfile(m,n) = h0(N*(m-1) + n);
    end
end
%--------------------------------------------------------------------------
function [cellsNew,mutant1New,mutant2New,mutant3New] = kUpdateCells(...
    cells,gProfile,cProfile,hProfile,phia,mutant1,mutant2,mutant3,a0,hN,...
    hT,pc,pa,k,nA)

[M,N] = size(cells);
cellsNew = zeros(M+1,N); cellsNew(1:M,:) = cells;
mutant1New = zeros(M+1,N); mutant1New(1:M,:) = mutant1;
mutant2New = zeros(M+1,N); mutant2New(1:M,:) = mutant2;
mutant3New = zeros(M+1,N); mutant3New(1:M,:) = mutant3;

P = kPermZmZn(M,N);

for X=P

    m = X(1); n = X(2);

    m1 = mutant1New(m,n);
    m3 = mutant3New(m,n);
    h = hProfile(m,n);
    a = phia(m,n);

    pdiv = (a - a0)/(1 - a0);

    if (m3 == 0)

        pdea = h/hN;

    else

        pdea = h/hT;

    end

    if (cells(m,n) == 0)

    elseif (cells(m,n) == -1)

        if (rand < pc)
            cellsNew(m,n) = 0;
        end

    elseif ((rand < pdea)|(a < a0)|((m1 == 0)&(m > 1)))

        cellsNew(m,n) = -1;
        mutant1New(m,n) = 0;
        mutant2New(m,n) = 0;
        mutant3New(m,n) = 0;

    else
        if (rand < pa)
            z = randsrc(1,1,[1,2,3]);

            if (z==1)

                mutant1New(m,n) = mod(1+mutant1New(m,n),2);

            elseif (z == 2)

                mutant2New(m,n) = mod(1+mutant2New(m,n),2);

            elseif (z == 3)

                mutant3New(m,n) = mod(1+mutant3New(m,n),2);
            end
        end

        if (rand < pdiv)

            C = zeros(1,4);

            if (cellsNew(m,kMod(n-1,N)) == 0)

                C(1) = cProfile(m,kMod(n-1,N));

            end

            if (cellsNew(m,kMod(n+1,N)) == 0)

                C(2) = cProfile(m,kMod(n+1,N));

            end

            if (cellsNew(m+1,n) == 0)

                if (m ~= M)

                    C(3) = cProfile(m+1,n);

                else

                    C(3) = cProfile(m,n);

                end
            end

            if (m ~= 1)

                if (cellsNew(m-1,n) == 0)

                    C(4) = cProfile(m-1,n);

                end
            end

            if (max(C) > 0)

                [temp,j] = max(C);

                if (j == 1)
                    m1 = m; n1 = kMod(n-1,N);
                elseif (j == 2)
                    m1 = m; n1 = kMod(n+1,N);
                elseif (j == 3)
                    m1 = m+1; n1 = n;
                elseif (j == 4)
                    m1 = m-1; n1 = n;
                end

                mut1 = mutant1New(m,n);
                mut2 = mutant2New(m,n);
                mut3 = mutant1New(m,n);

                h = hProfile(min(m1,M),n1);
                g = gProfile(min(m1,M),n1);
                c = cProfile(min(m1,M),n1);
                phig = (1 + (k-1)*mut2).*g;
                a = c + nA*(phig - c);

                if ((a < a0)|((mut1 == 0)&(m1 > 1))|(h > hT)|...
                        ((h > hN)&(mut3 == 0)))

                else
                    cellsNew(m1,n1) = 1;
                    mutant1New(m1,n1) = mutant1New(m,n);
                    mutant2New(m1,n1) = mutant2New(m,n);
                    mutant3New(m1,n1) = mutant3New(m,n);
                end

            end
        end
    end
end

if (cellsNew(M+1,:) == zeros(1,N))

    cellsNew = cellsNew(1:M,:);
    mutant1New = mutant1New(1:M,:);
    mutant2New = mutant2New(1:M,:);
    mutant3New = mutant3New(1:M,:);

end
%--------------------------------------------------------------------------
function gProfile = gApprox(deltag);

[M,N] = size(deltag);

K = sqrt(mean(mean(deltag)));
C1 = (1+exp(2*K*M))^(-1);
C2 = (tanh(K*M)+1)/2;

gProfile = zeros(M,N);

for m = 1:M
    gProfile(m,:) = C1*exp(K*m) + C2*exp(-K*m);
end
%--------------------------------------------------------------------------
function Xall = kUpdateAll(Xnew,Xall)

[M,N,z] = size(Xall);
[M1,temp] = size(Xnew);

if (isempty(Xall))

    Xall = Xnew;

elseif (M == M1)

    Xall(:,:,z+1) = Xnew;

else

    Xalltemp = zeros(M1,N,z+1);
    Xalltemp(1:M,:,1:z) = Xall;
    Xalltemp(:,:,z+1) = Xnew;
    Xall = Xalltemp;

end
%--------------------------------------------------------------------------