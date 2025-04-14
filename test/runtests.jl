using BinnedDistributionFit
using Test, FHist

@testset "ExtendPdf" begin
    d = ExtendPdf(x->x, (1,3))
    h = Hist1D(; binedges=1:3, bincounts=[2.0, 4.0], sumw2=[2.0, 4.0])
    NF = RooFitNLL_functor(d, h; num_integrator=BinnedDistributionFit.QuadGKIntegrator())
    @test NF(2.0) ≈ 1.6827899396467232
end
