.PHONY: run-daemon figures serve-docs

figures:
	poetry run brei figures

run-daemon:
	julia --project=. --startup-file=no -e 'using Revise; using DaemonMode; serve()'

serve-docs:
	julia --project=docs -e 'using LiveServer; servedocs()'

