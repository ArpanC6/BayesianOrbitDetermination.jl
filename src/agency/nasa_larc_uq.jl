# NASA LaRC Uncertainty Quantification (UQ) Challenge Problem Benchmark Solver

struct NASAUQChallengeResult
    epistemic_bounds::Tuple{Float64, Float64}
    aleatory_std::Float64
    reliability_index::Float64
end

"""
    solve_nasa_larc_uq_subproblem(r_samples::Matrix{Float64}, u_true::Vector{Float64}) -> NASAUQChallengeResult

Propagate mixed aleatory and epistemic uncertainties for NASA Langley Research Center UQ Challenge problem.
"""
function solve_nasa_larc_uq_subproblem(posterior_samples::Matrix{Float64}, u_true::Vector{Float64})
    pos_errors = [norm(posterior_samples[i, 1:3] - u_true[1:3]) for i in 1:size(posterior_samples, 1)]
    
    # Epistemic 95% confidence interval on position error
    q025 = quantile(pos_errors, 0.025)
    q975 = quantile(pos_errors, 0.975)
    
    aleatory_sigma = std(pos_errors)
    reliability = mean(pos_errors .< 0.1) # Fraction within 100m tolerance

    return NASAUQChallengeResult((q025, q975), aleatory_sigma, reliability)
end
