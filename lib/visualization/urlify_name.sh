#!/usr/bin/env bash

echo "$1" | iconv -t ASCII//TRANSLIT | sed -E "s/[ \/]/_/g; s/['\"]//g"
