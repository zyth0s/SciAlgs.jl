using StatsBase

"""Implements MCMC sampling for the hardcore model."""
function mcmc_hardcore(grid::BitArray{2}, trials::Int)
  (m, n) = size(grid)
  sites_occupied = zeros(Int, trials)
  total = 0
  for trial = 1:trials
    sites_occupied[trial] = total
    i = rand(1:m)
    j = rand(1:n)
    if grid[i, j]
      grid[i, j] = false
      total = total - 1
      sites_occupied[trial] = total
    else
      if i > 1 && grid[i - 1, j]
       continue
      end
      if i < m && grid[i + 1, j]
        continue
      end
      if j > 1 && grid[i, j - 1]
        continue
      end
      if j < n && grid[i, j + 1]
        continue
      end
      grid[i, j] = true
      total = total + 1
      sites_occupied[trial] = total
    end
  end
  return sites_occupied
end

(m, n, trials) = (50, 50, 1000000);
grid = falses(m, n);
sites_occupied = mcmc_hardcore(grid, trials);
describe(sites_occupied[10000:trials])
