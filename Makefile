.PHONY: run-daemon figures serve-docs docs test

figures:
	poetry run brei figures

run-daemon:
	julia --project=. --startup-file=no -e 'using Revise; using DaemonMode; serve()'

serve-docs:
	julia --project=docs -e 'using LiveServer; servedocs()'

docs:
	julia --project=docs -e 'include("docs/make.jl")'

test:
	julia --project=. -e 'using Pkg; Pkg.test()'

pluto:
	julia --project=workenv -e 'using Pluto; Pluto.run()'
