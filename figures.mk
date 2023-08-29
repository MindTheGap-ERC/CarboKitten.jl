# ~/~ begin <<docs/src/carbocat-ca.md#figures.mk>>[init]
.RECIPEPREFIX = >
.PHONY: all _all

fig := docs/fig

all: _all

# ~/~ begin <<docs/src/carbocat-ca.md#build>>[init]
targets += $(fig)/burgess2013-fig3.svg
targets += $(fig)/burgess2013-long-times.svg

docs/fig/burgess2013-fig3.svg: src/figures/ca.jl
> julia --project=. -e 'include("$<"); plot("$@")'

docs/fig/burgess2013-long-times.svg: src/figures/ca.jl
> julia --project=. -e 'include("$<"); plot_long_times("$@")'
# ~/~ end

_all: $(targets)
# ~/~ end