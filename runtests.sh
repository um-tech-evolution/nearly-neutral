#! /usr/bin/env bash

# Get the directory where the script is located. Should be the
# repository root directoy.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

usage() {
    echo "Usage: ./runtests.sh"
}

runtests() {
    echo "Running tests using random number seed 1:"
    cd "$DIR/src" && julia --color=yes --check-bounds=yes run.jl examples/ia_example1 1
    cd "$DIR/src" && julia --color=yes --check-bounds=yes run.jl examples/ia_example2 1
    cd "$DIR/src" && julia --color=yes --check-bounds=yes run.jl examples/ia_example3 1
    cd "$DIR/src" && julia --color=yes --check-bounds=yes run.jl examples/ia_example4 1
    cd "$DIR/src" && julia --color=yes --check-bounds=yes run.jl examples/ia_example5 1
    cd "$DIR/src" && julia --color=yes --check-bounds=yes run.jl examples/is_example0 1
    cd "$DIR/src" && julia --color=yes --check-bounds=yes run.jl examples/is_example1 1
    cd "$DIR/src" && julia --color=yes --check-bounds=yes run.jl examples/is_example2 1
}

runtests

exit 0

