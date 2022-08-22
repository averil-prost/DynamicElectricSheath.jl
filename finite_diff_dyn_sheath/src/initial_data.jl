export E0, fi_0, fe_0

"""
$(SIGNATURES)

Electric field
"""
E0(x) = 0.0

"""
$(SIGNATURES)

Initial mask
"""
mask(x) = 0.5 * (tanh((x + 0.1) ./ 0.1) - tanh((x - 0.1) ./ 0.1))

"""
$(SIGNATURES)

Ions distribution
"""
fi_0(x, v; μ = 0.5) = mask(x) * exp(-0.5 .* v^2) / sqrt(2π)

"""
$(SIGNATURES)

Electrons distribution
"""
fe_0(x, v; μ = 0.5) = sqrt(μ) * mask(x) * exp(-0.5 * μ * v^2) / sqrt(2π)
