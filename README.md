# BinnedDistributionFit

[![Build Status](https://github.com/Moelf/BinnedDistributionFit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Moelf/BinnedDistributionFit.jl/actions/workflows/CI.yml?query=branch%3Amain)

## [Work In Progress]
Everything is subject to change, including type names and type hierarchy.

## The Scope of this Package
We try to cover a core set of usages common in HEP, namely, fitting one or multiple distributions to binned data, using
Chi2  or Likelihood objective.

Conceptually, given user defined function ("shape"), we can treat it as if it's a PDF by numerically normalizing it at
every evaluating point.


## Math for One PDF fitting to Data: (i.e. `ExtendPdf`)

For main.py and main.jl (just the ExtendPdf)

if you change the observed bincounts from `[2.0, 4.0]` to `[4.0, 8.0]`, the NLL becomes `1.3655798792934464`

If on top of that, you also change `p0` from 2.0 to 4.0, you get NLL `-4.952186287425897`

![](https://gist.github.com/user-attachments/assets/87f7926f-9098-4b94-be74-98c6db2fe0d6)

## Math for sum of two PDFs fitting to data: (i.e. `SUM::` in RooFit)

![](https://gist.github.com/user-attachments/assets/b90aaf38-dbf8-49b4-b564-152f54b9cc38)
