#!/usr/bin/env bash

# converts a primer set name into a form usable in a URL
# utility script used in multiple places in the visualization scripts

echo "$1" | iconv -t ASCII//TRANSLIT | sed -E "s/[ \/]/_/g; s/['\"]//g"
