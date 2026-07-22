#!/bin/bash
set -o pipefail
export PATH="/root/miniconda3/envs/md_env/bin:$PATH"
source /root/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate md_env 2>/dev/null
cd /root/autodl-tmp/docking
MAX_RETRIES=10
for a in $(seq 1 $MAX_RETRIES); do
    echo "=== Attempt $a at $(date) ===" | tee -a runner.log
    PYTHONUNBUFFERED=1 python3 -u md_salt_bridge_cuda.py 2>&1 | tee -a runner.log
    rc=${PIPESTATUS[0]}
    [ $rc -eq 0 ] && echo "=== SUCCESS ===" | tee -a runner.log && exit 0
    echo "=== FAILED (rc=$rc) at $(date), retry ${a}/${MAX_RETRIES} ===" | tee -a runner.log
    nvidia-smi --query-gpu=timestamp,name,utilization.gpu,memory.used --format=csv,noheader >> runner.log
    sleep 60
done
echo "=== EXHAUSTED ===" | tee -a runner.log
exit 1
