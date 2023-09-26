#!/usr/bin/env bash

gvpack -u pipeline_documentation.dot | dot "$@"
