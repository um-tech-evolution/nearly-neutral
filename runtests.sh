#! /usr/bin/env bash

# Get the directory where the script is located. Should be the
# repository root directoy.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
    echo "Usage: ./runtests.sh"
}

runtests() {
    echo "Running tests:"
    #cd "$DIR/experiments" && julia --color=yes --check-bounds=yes -e "include(\"n_neutral.jl\");"
    #cd "$DIR/experiments" && julia --color=yes --check-bounds=yes -e "include(\"n_infsites.jl\");"
    #cd "$DIR/experiments" && ls -l
    #cd "$DIR/experiments" && julia --color=yes --check-bounds=yes n_neutral.jl
    cd "$DIR/experiments" && julia --color=yes --check-bounds=yes n_neutral.jl examples/nn_example1
    cd "$DIR/experiments" && julia --color=yes --check-bounds=yes n_neutral.jl examples/nn_example2
    cd "$DIR/experiments" && julia --color=yes --check-bounds=yes n_infsites.jl examples/in_example1
    cd "$DIR/experiments" && julia --color=yes --check-bounds=yes n_infsites.jl examples/in_example2
    #cd "$DIR/experiments" && julia --color=yes --check-bounds=yes n_fixed.jl examples/fix_example1
    #cd "$DIR/experiments" && julia --color=yes --check-bounds=yes n_fixed.jl examples/fix_example2
}

runtests

exit 0

