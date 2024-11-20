#!/bin/sh
#
# Entry point script for a docker image that simply runs the domain decomposition tool, e.g
# docker build . -t domain_decomp && docker run --rm -v $PWD:/io domain_decomp 8 -g testgrid.nc

. /opt/spack-environment/activate.sh

NP=$1; shift

case $NP in
    ''|*[!0-9]*)
        echo "Provide the number of domains before any options"
        exec decomp --help
        exit 15
        ;;
    *) mpirun --allow-run-as-root -n $NP --oversubscribe decomp "$@" ;;
esac
