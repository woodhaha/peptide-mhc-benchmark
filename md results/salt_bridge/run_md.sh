#!/bin/bash
source /root/miniconda3/etc/profile.d/conda.sh
conda activate md_env
cd /root/autodl-tmp/docking
PYTHONUNBUFFERED=1 nohup python3 -u md_salt_bridge_cuda.py > md_salt_bridge.log 2>&1 &
echo "PID=$! launched"
