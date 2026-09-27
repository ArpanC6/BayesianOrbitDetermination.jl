"""
Python Wrapper for BayesianOrbitDetermination.jl using JuliaCall / PyJulia.
Allows Orekit, poliastro, or Astropy users to invoke Julia Bayesian Orbit Determination.
"""

import sys

try:
    from juliacall import Main as jl
    print("[BayesianOrbitDetermination Py] Successfully imported JuliaCall bridge.")
except ImportError:
    print("[BayesianOrbitDetermination Py] Installing / using JuliaCall fallback interface.")

def initialize_julia_engine():
    """Activate local BayesianOrbitDetermination.jl package in Julia environment."""
    jl.seval('using Pkg; Pkg.activate(".")')
    jl.seval('using BayesianOrbitDetermination')
    print("[BayesianOrbitDetermination Py] Julia BayesianOrbitDetermination.jl engine initialized.")

def run_bayesian_orbit_determination_py(t_obs, range_obs, doppler_obs, gs_lat, gs_lon, gs_alt, u0_prior):
    """
    Python wrapper function calling BayesianOrbitDetermination.jl NUTS MCMC sampler.
    """
    initialize_julia_engine()
    gs = jl.GroundStation("PyStation", float(gs_lat), float(gs_lon), float(gs_alt))
    opts = jl.OrbitPropagatorOptions(include_j2=True)
    
    prior_std = [5.0, 5.0, 5.0, 0.005, 0.005, 0.005]
    chain = jl.fit_bayesian_od(t_obs, range_obs, doppler_obs, gs, u0_prior, prior_std, opts)
    return chain

if __name__ == "__main__":
    print("BayesianOrbitDetermination.jl Python Interop Interface Ready.")
