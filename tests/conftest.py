"""The tests import this repository's `report` package and run its scripts; the `ngsdose` library must be
installed (pip install -e ../NGS-DOSE) with its GRCh38 bundle reachable (NGSDOSE_RESOURCES, or the checkout)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
