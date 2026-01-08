.PHONY: run-daemon figures serve-docs docs test

figures:
	poetry run brei figures

run-daemon:
	julia --project=workenv --startup-file=no -e 'using Revise; using DaemonMode; serve()'

serve-docs:
	julia +1.12 --project=docs -e 'using Pkg; Pkg.instantiate(); using Revise; using LiveServer; servedocs()'

docs:
	julia +1.12 --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'

test:
	julia +1.10 --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

pluto:
	julia --project=workenv -e 'using Pluto; Pluto.run()'
