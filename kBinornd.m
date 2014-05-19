function r = kBinornd(n,p)
% Random number from the binomial distribution.


if nargin < 2
    error('stats:binornd:TooFewInputs','Requires at least two input arguments.'); 
end


% Handle the scalar params case efficiently
if isscalar(n) & isscalar(p) % scalar params
    if (0 <= p & p <= 1) & (0 <= n & round(n) == n)
        r = sum(rand([1,double(n)]) < p, 2);
    else
        r = repmat(NaN, 1);
    end

end
