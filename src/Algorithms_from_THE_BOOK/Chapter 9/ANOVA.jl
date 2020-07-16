"""Estimates ANOVA parameters with no missing data."""
function two_way_ANOVA(Y::Array{T, 3}) where T <: Real
  (i, j, k) = size(Y)
  S = sum(Y, dims = 3)
  (alpha, beta) = (sum(S, dims = 2), beta = sum(S, dims = 1))
  mu = sum(alpha) / (i * j * k)
  alpha = alpha / (j * k) .- mu
  beta = beta / (i * k) .- mu
  return (mu, alpha, beta)
end

"""Estimates ANOVA parameters with missing data by an MM
algorithm. The observed data appear in Y. The corresponding 
entries of W are 0 for missing data and 1 for observed data."""
function MM_ANOVA(Y::Array{T, 3}, W::Array{T, 3}) where T <: Real
  (i, j, k) = size(Y)
  (mu, alpha, beta) = (zero(T), zeros(T, i), zeros(T, j))
  P = zeros(size(Y)) # predicted values
  X = zeros(size(Y))
  old_rss = Inf
  for n = 1:100 # mm iteration loop
    rss = zero(T) # residual sum of squares
    for ii = 1:i
      for jj = 1:j
        for kk = 1:k
          P[ii, jj, kk] = mu + alpha[ii] + beta[jj]
          residual = Y[ii, jj, kk] - P[ii, jj, kk]
          rss = rss + W[ii, jj, kk] * residual^2
        end
      end
    end
    X = W .* Y + (one(T) .- W) .* P
    (mu, alpha, beta) = two_way_ANOVA(X) # MM update
    if abs(old_rss - rss) < 1e-8 # check for convergence
      return (mu, alpha, beta)
    end
    old_rss = rss
  end
  return (mu, alpha, beta)
end

Y = zeros(3, 2, 4);
Y[:, :, 1] = [1 2; 3 4; 5 6];
Y[:, :, 2] = -Y[:, :, 1];
Y[:, :, 3] = Y[:, :, 1] / 2; 
Y[:, :, 4] = -Y[:,:,1] / 2; 
W = ones(size(Y));
W[1, 1, 1] = 0.0;
W[3, 2, 4] = 0.0;
(mu, alpha, beta) = MM_ANOVA(Y, W);
println("mu = ",mu);
println("alpha = ",alpha);
println("beta = ",beta);



  