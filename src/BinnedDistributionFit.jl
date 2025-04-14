module BinnedDistributionFit

export ExtendPdf, RooFitNLL, RooFitNLL_functor, pdf

import NumericalDistributions: NumericallyIntegrable, pdf, cdf
using FHist: Hist1D, binedges, bincenters, bincounts
using QuadGK: quadgk    
using Trapz: trapz

"""
    struct ExtendPdf{F, S}
        func::F
        support::S
    end

Type representing an extended pdf where normalization does not need to equal to 1.
`func` can be any function that can be evaluated within `support`, when calculating Maximum Likelihood, `func` will be numerically normalized over `support`.

## Examples
```julia
julia> f(x, ps) = x^ps[1] + x*ps[2]

julia> d = BinnedDistributionFit.ExtendPdf(f, (1, 2))

julia> BinnedDistributionFit.pdf(d, 2.0, [1,2])
6.0

julia> BinnedDistributionFit.pdf(d, 2.0, [1,2]; num_integral=2.0)
3.0
```
"""
struct ExtendPdf{F, S}
    func::F
    support::S
end

include("./num_integrator.jl")
include("./RooFitNLL.jll")

end
