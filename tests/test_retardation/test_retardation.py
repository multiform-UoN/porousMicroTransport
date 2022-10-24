import pytest

import subprocess
from pathlib import Path

import numpy as np

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

DIR = Path(__file__).parent

@pytest.fixture(scope="module")
def retardation_case():
    subprocess.run(["./clean"],cwd=DIR)
    subprocess.run(["./run"], check=True, cwd=DIR)

def test_retardation(retardation_case):
    for d in DIR.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue

        if t==0:
            continue

        A = np.array(ParsedParameterFile(DIR / d.name / "A")["internalField"].value())
        B = np.array(ParsedParameterFile(DIR / d.name / "B")["internalField"].value())

        # Check total mass
        assert np.sum(B) == pytest.approx(np.sum(A)/2, rel=1e-3)

        # Check distributions
        assert B[:len(B)//2] == pytest.approx(A[::2], abs=1e-1)
