function [xhat, P] = condition(jointdist, y, condind)
% CONDITION approximate the conditional mean and covariance for a joint distibution
%
%  [xhat, P] = condition(jointdist, y, condInd)
%
%  jointdist   the distribution for which to approximate the conditional
%              mean and covariance.
%  y           value to condition on.
%  condind     an index vector of the components that should be conditioned
%              on.
%  xhat        the approximated conditional mean
%  P           the approximated conditional covariance

% Copyright Gustaf Hendeby
%$ Revision: 28-Oct-2019 $

  if ~isa(jointdist, 'pdfclass')
    error('The first argument must be a distribution');
  end
  mu = mean(jointdist);
  Sigma = cov(jointdist);

  if islogical(condind)
    if sum(condind) >= numel
      error('The third argument must match a subset of indices in jointdist');
    end
    I = condind;
  else
    if numel(condind) >= numel(mu)
      error('The third argument must match a subset of indices in jointdist');
    end
    I =  ismember(1:numel(mu), condind);
  end

  xhat = mu(~I);
  yhat = mu(I);
  Pxx = Sigma(~I, ~I);
  Pyy = Sigma(I, I);
  Pxy = Sigma(~I, I);

  xhat = xhat + (Pxy/Pyy)*(y-yhat);
  P = Pxx - (Pxy/Pyy)*Pxy';
end
