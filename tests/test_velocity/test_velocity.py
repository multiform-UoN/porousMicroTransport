import pytest

import subprocess
from pathlib import Path

import numpy as np

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

DIR = Path(__file__).parent

@pytest.fixture(scope="module")
def velocity_case():
    subprocess.run(["./clean"],cwd=DIR)
    subprocess.run(["./run"], check=True, cwd=DIR)

def test_infiltration(velocity_case):
    theta0 = np.array(ParsedParameterFile(DIR / "0" / "theta")["internalField"].value())
    U = np.array(ParsedParameterFile(DIR / "0" / "U")["boundaryField"]["left"]["value"].value())[0]
    h = 30e-3/5000

    for d in DIR.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue

        theta = np.array(ParsedParameterFile(DIR / d.name / "theta")["internalField"].value())

        assert np.sum(theta - theta0)*h == pytest.approx(U*t)
