export trial_result

@doc """ type trial_result
  Julia type that stores both the parameters and the results of a trial.
"""
mutable struct trial_result
  nn_simtype::Int64
  n::Int64    # sample size, must be <= N
  N::Int64    # popsize
  N_mu::Float64  # N*mu, population mutation rate.  theta = 2*N_mu
  ngens::Int64
  burn_in::Float64  # The number of generations of burn-in is N/mu + 50
  dfe::Function
  dfe_str::AbstractString
  use_poplist::Bool
  mean_fitness::Float64
  expected_richness::Float64  # sum_{i=0}^{n-1} theta/(theta+i) where theta = 2*N_mu.
  average_richness::Float64   # Average number of traits in populations
  stderr_richness::Float64
  expected_w_homoz::Float64
  w_homoz::Float64
  stderr_w_homoz::Float64
  IQV::Float64
  stderr_IQV::Float64
  entropy::Float64
  stderr_entropy::Float64
end

# Constructor that sets the parameters
function trial_result( nn_simtype::Int64, n::Int64, N::Int64, N_mu::Float64, ngens::Int64,
    burn_in::Float64=2.0, dfe::Function=dfe_neutral,
    dfe_str::AbstractString="neutral"; use_poplist=false )
  tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, dfe, dfe_str, use_poplist, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
  add_expected_richness( tr )
  tr
end

