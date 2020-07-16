using DecisionTree, LinearAlgebra, SparseArrays

"""Performs matrix completion discriminant analysis. M is the sparse
predictor matrix. Missing entries of class are coded as 0."""
function mcda(M::AbstractMatrix{T}, class::Vector{Int}, 
  classes::Int, r::Int, epsilon::T) where T <: Real
#
  T = eltype(epsilon)
  (cases, features) = size(M)
  vertex = zeros(T, classes, classes - 1)
  (b, c, d) = (classes - 1, sqrt(classes), sqrt(classes - 1))
  vertex[1, :] .= one(T) / d
  for i = 1:(classes - 1) # create the simplex vertices
    vertex[i + 1, :] .= -(one(T) + c) / (b * d)
    vertex[i + 1, i] = vertex[i + 1, i] + c / d
  end
  Z = zeros(T, cases, classes - 1) # discriminant columns
  for case = 1:cases 
    if class[case] != 0
      Z[case, :] .= vertex[class[case], :]
    end
  end
  Y = [M sparse(Z)]
  (maxiters, tol) = (1000, 1e-6) # carry out MCDA
  R = similar(Y)
  U = randn(cases, r)
  V = randn(r, features + classes - 1)
  old_loss = 1e20
  for iteration = 1:maxiters # use alternating least squares and MM
    for i in eachindex(Y) # compute current residuals
      (j, k) = i.I
      R[i] = Y[i] - dot(U[j, :], V[:, k])
    end
    S = (U' * U + epsilon * I) \ U' # update V
    V = S * R + (S * U) * V
    for i in eachindex(Y) # compute current residuals
      (j, k) = i.I
      R[i] = Y[i] - dot(U[j, :], V[:, k])
    end
    T = (V * V' + epsilon * I) \ V # update U
    U = R * T' + U * (V * T')
    loss = norm(R)^2 # check for convergence
    if abs(old_loss - loss) < (old_loss + 1.0) * tol
      break
    end
    old_loss = copy(loss)
  end
  assigned = copy(class) # classify the unclassified cases
  d = zeros(classes)
  for case = 1:cases
    if class[case] == 0
      predicted = vec(U[case, :]' * V[:, features + 1:end])
      for j = 1:classes
        d[j] = norm(predicted  - vec(vertex[j, :]))
      end
      assigned[case] = argmin(d) 
    end
  end
  return assigned
end

"""Prepares Fisher's Iris data for analysis."""
function analyze_iris_data(epsilon::T, missing_rate::T,
  classes::Int) where T <: Real
#
  (features, labels) = load_data("iris")  # data input Array{Any}
  M = float.(features)
  label = string.(labels)
  cases = size(M, 1)
  class = zeros(Int, cases)
  for i = 1:cases
    if label[i] == "Iris-setosa"
      class[i] = 1
    elseif label[i] == "Iris-versicolor"
      class[i] = 2
    elseif label[i] == "Iris-virginica"
      class[i] = 3
    else
      println("failed label")
    end
  end
  for i in eachindex(M) # randomly delete features
    if rand() < missing_rate
      M[i] = 0.0
    end
  end
  M = sparse(M)
  deletedclass = copy(class)
  for i = 1:cases # randomly delete classes
    if rand() < missing_rate
      deletedclass[i] = 0
    end
  end
  for r = 1:6 # run MCDA
    imputedclass = mcda(M, deletedclass, classes, r, epsilon);
    errors = count(class - imputedclass != 0);
    println("rank = ",r," classification_errors = ", errors)
  end
  return nothing
end

(classes, epsilon, missing_rate) = (3, 1e-5, 0.25);
analyze_iris_data(epsilon, missing_rate, classes) 
