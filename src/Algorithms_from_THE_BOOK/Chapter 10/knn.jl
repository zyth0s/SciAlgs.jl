using Distances, StatsBase

"""Performs k nearest neighbor classification with training data
Y. The classes should be numbered 1, 2,...""" 
function knn(X::Matrix{T}, Y::Matrix{T}, class::Vector{Int}, 
  k::Int) where T <: Real
#
  testing = size(X, 2)
  predicted_class = zeros(Int, testing)
  distance = pairwise(Euclidean(), Y, X, dims = 2)
  for i = 1:testing # find k nearest neighbors
    perm = partialsortperm(distance[:, i], 1:k)
    predicted_class[i] = mode(class[perm]) # most common class
  end
  return predicted_class
end

(training, testing, features) = (100, 10, 30);
X = randn(features, testing);
Y = randn(features, training);
(k, classes) = (3, 2);
class = rand(1:classes, training);
predicted_class = knn(X, Y, class, k)
