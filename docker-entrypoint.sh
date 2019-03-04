#!/bin/bash
umask 0002
. /opt/conda/etc/profile.d/conda.sh
conda activate betaduck
betaduck "${@}"
