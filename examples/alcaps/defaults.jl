# ~/~ begin <<docs/src/model-alcap.md#examples/alcaps/defaults.jl>>[init]
#| requires: src/Model/ALCAPS.jl
#| creates: data/alcaps_default.h5

using CarboKitten.Model.ALCAPS

ALCAPS.main(ALCAPS.Input(), "data/alcaps_default.h5")
# ~/~ end