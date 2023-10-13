#!/usr/bin/env bash

# converts a primer set name into a form usable in a URL

echo "$1" | iconv -t ASCII//TRANSLIT | sed -E "s/[ \/]/_/g; s/['\"]//g"
