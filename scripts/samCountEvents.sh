#!/bin/bash
# Evaluate N(gen events)
if [[ "$1" == "" ]]; then
    exit
fi

samweb list-files --summary "dh.dataset=${1} and availability:anylocation" | while read line
do
  read -a strarr <<< "$line"
  if [[ "${strarr[0]}" =~ "Event" ]]; then
    echo "${strarr[2]}"
  fi
done
