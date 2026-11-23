#!/usr/bin/env python3
"""Run MAFtools from a checkout without installing it."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from maftools.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
