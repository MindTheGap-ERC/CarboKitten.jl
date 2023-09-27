.PHONY: run-daemon

figures:
	make -f .entangled/build/Makefile

run-daemon:
	julia --project=. --startup-file=no -e 'using DaemonMode; serve()'

serve-docs:
	julia --project=docs -e 'using LiveServer; servedocs()'

