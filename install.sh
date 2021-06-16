#!/usr/bin/env bash -ue

function link_to {
    local dst=$1
    local pwd=$(readlink -f $PWD)
    for f in $(ls *.py); do
        command chmod +x $pwd/$f
        command ln -snf $pwd/$f $local_bin/$(echo $f | sed -e "s@\.py@@")
    done
}

function main {
    local local_bin=$HOME/.local/bin
    if ! [[ -e $local_bin ]]; then
        command echo "${local_bin} is not exists, auto make it."
        command mkdir -p $local_bin
    fi

    link_to local_bin
}

main
