.PHONY: run-daemon figures serve-docs

figures:
	make -f .entangled/build/Makefile

run-daemon:
	julia --project=. --startup-file=no -e 'using Revise; using DaemonMode; serve()'

serve-docs:
	julia --project=docs -e 'using LiveServer; servedocs()'

