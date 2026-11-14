# Zonal Earth gravity harmonics (J2, J3, J4, J5).
#
# The acceleration is the gradient of
#     U = (mu / r) * sum_n J_n (R_e / r)^n P_n(z / r),
# with the sign convention a = -grad U. Each term is checked against
# j2_perturbation for n = 2 and against a finite-difference gradient of U in
# the test suite.

const J3_EARTH = -2.5324105e-6
const J4_EARTH = -1.6199216e-6
const J5_EARTH = -2.2788820e-7

const ZONAL_J = (J2_EARTH, J3_EARTH, J4_EARTH, J5_EARTH)

"""
    zonal_harmonic_perturbation(r; mu = MU_EARTH, R_e = R_EARTH,
                                n_max = 4, J = ZONAL_J)

Acceleration from the zonal harmonics J2 through `n_max` (default 4) in the
ECI frame. `J[n-1]` is the coefficient for degree `n`.
"""
function zonal_harmonic_perturbation(r::SVector{3, T};
                                     mu::Float64 = MU_EARTH,
                                     R_e::Float64 = R_EARTH,
                                     n_max::Int = 4,
                                     J::NTuple{4, Float64} = ZONAL_J) where T
    2 <= n_max <= 1 + length(J) || throw(ArgumentError("n_max must be in 2:$(1 + length(J))"))
    x, y, z = r[1], r[2], r[3]
    rn = norm(r)
    u = z / rn

    a = zero(SVector{3, T})
    P_prev2 = one(T)   # P_0
    P_prev1 = u        # P_1
    dP_prev2 = zero(T) # P_0'
    dP_prev1 = one(T)  # P_1'

    for n in 2:n_max
        Pn, dPn = if n == 2
            ((3u^2 - 1) / 2, 3u)
        else
            (((2n - 1) * u * P_prev1 - (n - 1) * P_prev2) / n,
             (2n - 1) * P_prev1 + dP_prev2)
        end

        coeff = mu / rn^2 * (R_e / rn)^n * J[n - 1]
        common = (n + 1) * Pn + u * dPn
        a += SVector{3, T}(
            coeff * (x / rn) * common,
            coeff * (y / rn) * common,
            coeff * ((n + 1) * u * Pn - (1 - u^2) * dPn),
        )

        P_prev2, P_prev1 = P_prev1, Pn
        dP_prev2, dP_prev1 = dP_prev1, dPn
    end
    return a
end
