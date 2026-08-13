"""
config.py - the only file you need to edit to run Stages 2 and 2b yourself.

Stages 2 and 2b download the satellite covariate stack from Google Earth
Engine. Everything Earth Engine needs to know about *you* lives in this file.
Everything the *manuscript* needs to stay reproducible lives further down, under
FROZEN, and must not be changed if you are reproducing published results.

Stages 1 and 3-9 do not use this file at all. If you obtained the pre-computed
covariate stack from the authors, you can ignore this file entirely.

Quick start:
    1. Set GEE_PROJECT below to your own Google Cloud project ID.
    2. Run:  earthengine authenticate
    3. Run:  python notebooks/stage2_covariate_stack.py
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent


# ===========================================================================
#  EDIT THIS
# ===========================================================================

# Your Google Cloud project ID, as a string.
#
# Where to find it: https://console.cloud.google.com/ - the project selector at
# the top of the page. It looks like "my-project-464812" or "ee-yourname".
# It is the *ID*, not the display name; they are often different.
#
# The project must have the Earth Engine API enabled, and your Google account
# must be registered for Earth Engine at https://earthengine.google.com/signup/
#
# Leave as None only if you are using a service-account key (see below), which
# carries its own project ID.
GEE_PROJECT: str | None = None      # e.g. "ee-yourname" or "my-project-464812"


# How to log in. Two options:
#
#   "auto"    - (default) use a service-account key if one is present at the
#               path below, otherwise fall back to the browser login.
#   "browser" - always use the browser login. Run `earthengine authenticate`
#               once, and Earth Engine remembers you on this machine.
#   "service" - always use the service-account key below. Used by the original
#               authors for unattended runs.
GEE_AUTH_MODE: str = "auto"


# Only relevant when GEE_AUTH_MODE is "service" or "auto".
#
# SECURITY: this is a credential, not data. It is git-ignored. Never commit it,
# never email it, never put it in a Zenodo archive. Anyone holding this file can
# spend money on your Google Cloud account.
GEE_SERVICE_ACCOUNT_KEY: Path = ROOT / ".secrets" / "gee_service_account.json"


# ===========================================================================
#  FROZEN - do not change these if you are reproducing the published results
# ===========================================================================
#
# These values are stated in the manuscript and are what make a rerun
# comparable to the published figures. They are recorded here for reference;
# the authoritative copies live at the top of each stage script.
#
#   Target CRS ................. EPSG:32636 (UTM 36N)
#   Grid resolution ............ 30 m
#   Sentinel-2 cropland mask ... date-locked to March 2025, to match the
#                                993-point ground survey window
#   Random seed ................ RANDOM_STATE = 42 in every stage
#   Monte Carlo draws .......... 200 (stage 5), 100 per scenario (stage 6)
#
# Changing any of the above produces a valid analysis of your own, but it is no
# longer a reproduction of El-Naggar et al. (2026). Say so if you publish it.
#
# A further caveat, which is not a bug and cannot be engineered away: Google
# periodically reprocesses the Sentinel archives. A stack rebuilt today will be
# very close to, but not bit-identical with, the stack used for the manuscript.
# To verify the exact published numbers, use the archived covariate stack from
# the companion Zenodo dataset rather than rebuilding it.


# ===========================================================================
#  Machinery - you should not need to read past this line
# ===========================================================================

def init_ee(verbose: bool = True):
    """Initialise the Earth Engine API according to the settings above.

    Raises RuntimeError with an actionable message rather than letting Earth
    Engine's own errors surface, because they are famously unhelpful.
    """
    try:
        import ee
    except ImportError as exc:                                # pragma: no cover
        raise RuntimeError(
            "The earthengine-api package is not installed.\n"
            "    Fix:  conda activate bani-suef-suitability\n"
            "    If that does not help:  pip install earthengine-api"
        ) from exc

    mode = GEE_AUTH_MODE.lower().strip()
    if mode not in {"auto", "browser", "service"}:
        raise RuntimeError(
            f"config.py: GEE_AUTH_MODE is {GEE_AUTH_MODE!r}, which is not one "
            "of 'auto', 'browser' or 'service'."
        )

    have_key = GEE_SERVICE_ACCOUNT_KEY.is_file()

    if mode == "service" and not have_key:
        raise RuntimeError(
            "config.py: GEE_AUTH_MODE is 'service' but no key file was found "
            f"at\n    {GEE_SERVICE_ACCOUNT_KEY}\n"
            "Either place the key there, or set GEE_AUTH_MODE = \"browser\" "
            "and run:  earthengine authenticate"
        )

    # --- service-account route --------------------------------------------
    if have_key and mode in {"auto", "service"}:
        import json
        try:
            info = json.loads(GEE_SERVICE_ACCOUNT_KEY.read_text())
            email, project = info["client_email"], info["project_id"]
        except (ValueError, KeyError) as exc:
            raise RuntimeError(
                f"The key file at {GEE_SERVICE_ACCOUNT_KEY} is not a valid "
                "Google service-account JSON key (expected 'client_email' and "
                "'project_id' fields)."
            ) from exc

        ee.Initialize(
            credentials=ee.ServiceAccountCredentials(email, str(GEE_SERVICE_ACCOUNT_KEY)),
            project=project,
        )
        if verbose:
            print(f"Earth Engine initialised via service account, project: {project}")
        return ee

    # --- browser route -----------------------------------------------------
    if not GEE_PROJECT:
        raise RuntimeError(
            "No Earth Engine credentials configured.\n\n"
            "  Open config.py in the project folder and set GEE_PROJECT to your\n"
            "  own Google Cloud project ID, for example:\n\n"
            '      GEE_PROJECT = "ee-yourname"\n\n'
            "  Then run this once, and follow the browser prompt:\n\n"
            "      earthengine authenticate\n\n"
            "  You need an Earth Engine account: https://earthengine.google.com/signup/\n"
            "  Approval usually takes 1-2 days.\n\n"
            "  To skip Stages 2 and 2b entirely, obtain the pre-computed stack\n"
            "  from the authors - see SETUP.md, Step 0, Bundle B."
        )

    try:
        ee.Initialize(project=GEE_PROJECT)
    except Exception as exc:
        raise RuntimeError(
            f"Earth Engine refused to initialise with project {GEE_PROJECT!r}.\n\n"
            "  Most common causes, in order:\n"
            "    1. You have not logged in yet.   Fix:  earthengine authenticate\n"
            "    2. GEE_PROJECT is the project *name*, not its *ID*. Check the\n"
            "       project selector at https://console.cloud.google.com/\n"
            "    3. The Earth Engine API is not enabled on that project.\n"
            "    4. Your Google account is not registered for Earth Engine:\n"
            "       https://earthengine.google.com/signup/\n\n"
            f"  Earth Engine said: {exc}"
        ) from exc

    if verbose:
        print(f"Earth Engine initialised via browser login, project: {GEE_PROJECT}")
    return ee


if __name__ == "__main__":
    # Lets a user test their credentials without running a 2-hour download:
    #     python config.py
    import sys

    try:
        ee = init_ee()
        n = ee.Number(1).add(1).getInfo()
        print("Round-trip test:", "OK" if n == 2 else f"unexpected result {n}")
        print("\nYou are ready to run stage 2:")
        print("    cd notebooks")
        print("    python stage2_covariate_stack.py")
        print("    python stage2b_download_stack.py")
    except RuntimeError as exc:
        print(f"\n{exc}\n", file=sys.stderr)
        sys.exit(1)
