"""Every shell script in the pipeline and the pilot parses."""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_shell_scripts_parse():
    scripts = sorted(list((ROOT / "pipeline").glob("*.sh")) + list((ROOT / "pilot").glob("*.sh")) + [ROOT / "regenerate.sh"])
    assert len(scripts) >= 10
    for sh in scripts:
        subprocess.run(["bash", "-n", str(sh)], check=True)
