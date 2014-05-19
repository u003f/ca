function r = kGeornd(p)
% Random number from the geometric distribution.

if nargin < 1
    error('stats:geornd:TooFewInputs','Requires at least one input argument.');
end



% Return NaN for elements corresponding to illegal parameter values.
p(p <= 0 | p > 1) = NaN;

% log(1-u) and log(1-p) "really" should both be negative, and
% the abs() here keeps the correct sign for their ratio when
% roundoff makes one of them exactly zero.
r = ceil(abs(log(rand(1)) ./ log(1 - p)) - 1); % == geoinv(u,p)

% Force a zero when p==1, instead of -1.
r(r < 0) = 0;
