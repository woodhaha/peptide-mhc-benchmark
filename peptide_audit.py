#!/usr/bin/env python3
"""Peptide-MHC manuscript consistency audit.

Thin wrapper over the shared config-driven engine. Edit `audit_config.json`
(not this file) when a number/name changes — those are the drift tripwires.
"""
import subprocess
import sys
from pathlib import Path

ENGINE = Path.home() / ".claude" / "skills" / "audit" / "audit_consistency.py"
CONFIG = Path(__file__).parent / "audit_config.json"

sys.exit(subprocess.call([sys.executable, str(ENGINE), str(CONFIG)]))
