#!/usr/bin/env python3
"""Fix the syntax error and re-run."""
import sys
path = "/root/autodl-tmp/docking/_extract_centroids.py"
src = open(path).read()
src = src.replace('print(f"\n{=*55}")', "print('='*55)")
src = src.replace('print(f"{=*55}")', "print('='*55)")
open(path, "w").write(src)
print("FIXED")
