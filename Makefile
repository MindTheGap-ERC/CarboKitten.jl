.PHONY: run-daemon figures serve-docs docs

figures:
	poetry run brei figures

run-daemon:
	julia --project=. --startup-file=no -e 'using Revise; using DaemonMode; serve()'

serve-docs:
	julia --project=docs -e 'using LiveServer; servedocs()'

docs:
	julia --project=docs -e 'include("docs/make.jl")'
