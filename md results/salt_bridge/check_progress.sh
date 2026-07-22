#!/bin/bash
# Check docking MD progress
LOG=/root/autodl-tmp/docking/md_salt_bridge.log
CSV=/root/autodl-tmp/docking/md_log.csv
echo "=== LOG (00:59:08) ==="
tail -5 "$LOG" 2>/dev/null
echo ""
if [ -f "$CSV" ]; then
    echo "=== CSV (last 3 lines) ==="
    tail -3 "$CSV"
    lines=$(wc -l < "$CSV")
    echo "Total lines: $lines"
fi
echo ""
echo "=== PROCESS ==="
ps -p $(pgrep -f md_salt_bridge_cuda) -o pid,pcpu,cputime,etime --no-headers 2>/dev/null || echo "Not running"
echo "=== GPU MEM ==="
nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader 2>/dev/null | grep -E "$(pgrep -f md_salt_bridge_cuda | head -1)" || echo "No GPU"
