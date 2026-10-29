"""Path resolution for the revision stages.

These scripts are distributed both inside the manuscript revision folder and
inside the code repository, so the repository root is located by walking up
from the script until a directory containing outputs/stage5 is found.
"""
from __future__ import annotations
import pathlib


def repo_root(start: pathlib.Path) -> pathlib.Path:
    here = pathlib.Path(start).resolve()
    for cand in [here, *here.parents]:
        if (cand / "outputs" / "stage5").is_dir():
            return cand
        for sub in cand.glob("*/outputs/stage5"):
            return sub.parents[1]
    raise FileNotFoundError(
        "repository root not found: expected a directory containing outputs/stage5")


REPO = repo_root(__file__)
OUT = REPO / "outputs" / "stage_rev"
WAPOR = REPO / "outputs" / "wapor"
OUT.mkdir(parents=True, exist_ok=True)
