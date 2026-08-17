#!/usr/bin/env python3
"""Backward-compatible entry point for the recovery-aware MCFR analysis.

The earlier hard-coded leakage-efficiency model has been superseded by
``paper2_mcfr_recovery_aware.py``.  This filename is retained so existing
reproduction commands continue to work.
"""
from paper2_mcfr_recovery_aware import main


if __name__ == "__main__":
    main()
