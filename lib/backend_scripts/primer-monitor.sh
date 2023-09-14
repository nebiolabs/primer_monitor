#!/usr/bin/env bash

SGE_QMASTER_PORT="$HEAD_SGE_QMASTER_PORT" SGE_ROOT="$HEAD_SGE_ROOT" "$QSUB_PATH" -S /bin/bash \
-v "NXF_JAVA_HOME=$NXF_JAVA_HOME" -v "BACKEND_INSTALL_PATH=$BACKEND_INSTALL_PATH" \
-o "$BACKEND_SCRATCH_PATH/download_qsub.log" -e "$BACKEND_SCRATCH_PATH/download_qsub.err" \
"$BACKEND_INSTALL_PATH}/lib/backend_scripts/primer-monitor-run.sh"