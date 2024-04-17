#!/bin/bash

scatterArg="$1"
shift
commandArgs="$@"

command="$(echo ${commandArgs} | sed "s|@@|${scatterArg}|g")"
${command}
