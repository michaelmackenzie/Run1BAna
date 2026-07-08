#!/bin/bash
# Return a dataset fileset
if [[ "$1" == "" ]]; then
    exit
fi

samweb list-files "dh.dataset=${1} and availability:anylocation"
