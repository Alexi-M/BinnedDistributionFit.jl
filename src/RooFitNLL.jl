"""
    RooFitNLL(normalized_predictions, observations; normalization)

Computes the negative log-likelihood for a binned distribution fit.
- `normalized_predictions`: A vector of normalized prediction at each bin-center.
- `observations`: A vector of observed values at each bin-center.
- `normalization`: A scalar value for normalization. Conceptually equal to the pdf generateing `predictions` integrated over the domain.

"""
function RooFitNLL(normalized_predictions, observations; normalization)
    per_bin_term = sum(@. observations * log(normalized_predictions))
    extend_term = sum(observations)*log(normalization) - normalization

    return -per_bin_term - extend_term
end

"""
    RooFitNLL_functor(d::ExtendPdf, data_hist::Hist1D; num_integrator=SimpleSumIntegrator())

Given an pdf and a data histogram, return a callable function that evaluates to the negative log likelihood with respect to the data histogram within the `d.support` domain.

The callable function should be called with N+1 parameter, the +1 being the first argument represening the overall # of events. For example, if user-defined function is `f(x, p) = x*p[1] + x^2*p[2]`, then the functor should be called with `[norm, p1, p2]`.

## Examples

```julia
julia> f(x,p) = x + p[1] + p[2]*x

julia> d = ExtendPdf(f, (1,3))

julia> BinnedDistributionFit.scalar_eval(d, 1, [2,3])
6

julia> NF = RooFitNLL_functor(d, h);

# the 10. is overall normalization
# [2.0, 3.0] are passed as `p` to `f(x,p)`
julia> NF([10., 2.0, 3.0])
0.060373400847997694
```
"""
function RooFitNLL_functor(d::ExtendPdf, data_hist::Hist1D; num_integrator=SimpleSumIntegrator())
    bes, bcs = binedges(data_hist), bincenters(data_hist)
    observations = bincounts(data_hist)
    @assert extrema(bes) == extrema(d.support) "The support of the distribution must match the bin edges of the histogram."

    function (norm_and_ps; kw...)
        norm, ps... = norm_and_ps
        predictions = vector_eval(d, bcs, ps; kw...)

        oneD_func(x) = scalar_eval(d, x, ps; kw...)
        numerical_int = _integrate(oneD_func, data_hist, num_integrator)
        normalized_predictions = predictions ./ numerical_int
        return RooFitNLL(normalized_predictions, observations; normalization=norm)
    end
end

"""
    RooFitNLL_functor(sumpdfs::SumOfPdfs, data_hist::Hist1D; num_integrator=SimpleSumIntegrator())

Given an `SumOfPdfs` and a data histogram, return a callable function that evaluates to the negative log likelihood with
respect to the data histogram within the `d.support` domain.

The callable function should be called with a vector of vector as input, specifically `[norms, ps1, ps2, ...]`, where `norms` is the normalization for each of the pdfs in order, and `ps1, ps2, ...` are the parameters for each pdf function, also in order.

## Examples

```julia
julia> f(x,ps) = x + ps[1] + ps[2]*x
julia> g(x,ps) = ps[1] + ps[3]

julia> d1 = ExtendPdf(f, (1,3))

julia> d2 = ExtendPdf(g, (1,3))

julia> sd = d1+d2

julia> h = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
edges: [1.0, 2.0, 3.0]
bin counts: [2.0, 4.0]
total count: 6.0

julia> NF = RooFitNLL_functor(sd, h; num_integrator=BinnedDistributionFit.QuadGKIntegrator())

julia> NF([[2.0, 3.0], [6.0, 7.0], [1.0, 5.0, 9.0]])
-0.627546320921633
```
"""
function RooFitNLL_functor(sumpdfs::SumOfPdfs, data_hist::Hist1D; num_integrator=SimpleSumIntegrator())
    bes, bcs = binedges(data_hist), bincenters(data_hist)
    observations = bincounts(data_hist)
    @assert extrema(bes) == extrema(sumpdfs.support) "The support of the distribution must match the bin edges of the histogram."

    Npdfs = length(sumpdfs.pdfs)
    function (norms_and_vps; kw...)
        # both are vectors
        norms, vps... = norms_and_vps
        overall_norm = sum(norms)
        fractions_of_pdfs = norms ./ overall_norm
        if length(norms) != Npdfs
            throw(ArgumentError("Expected $Npdfs normalizations, got $(length(norms))"))
        end

        integrals_of_pdfs = map(sumpdfs.pdfs, vps) do d, ps
            oneD_func(x) = scalar_eval(d, x, ps; kw...)
            _integrate(oneD_func, data_hist, num_integrator)
        end
        predictions_of_pdfs =  map(sumpdfs.pdfs, vps) do d, ps
            vector_eval(d, bcs, ps; kw...)
        end

        normalized_predictions = sum(fractions_of_pdfs .* predictions_of_pdfs ./ integrals_of_pdfs)
        return RooFitNLL(normalized_predictions, observations; normalization=overall_norm)
    end
end
