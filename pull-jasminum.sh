#!/bin/bash

REMOTE="scosta@10.101.72.16"
REMOTE_DIR="~/shape-optimizer"
LOCAL_DIR="."
PORT="27182"
KEY="~/.ssh/jasminum"

if [ -z "$1" ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

scp -r -P "${PORT}" -i "${KEY}" "${REMOTE}:${REMOTE_DIR}/$1" "${LOCAL_DIR}"