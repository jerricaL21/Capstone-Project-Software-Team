"""
install_pymol.py

Ensures the PyMOL Python API (pymol2 / pymol) is available in the SAME
Python interpreter that is currently running this script.

Why this matters:
    A very common failure mode is having PyMOL installed in one environment
    (e.g. a conda env, or a standalone PyMOL.app) while your actual app runs
    under a different interpreter (e.g. a system Python or venv). This script
    always targets `sys.executable`, so it installs into whatever interpreter
    is running it — no hardcoded paths, works anywhere.

Usage:
    python install_pymol.py            # install if missing, then verify
    python install_pymol.py --force    # reinstall/upgrade even if present
    python install_pymol.py --check    # only check, don't install

You can also import and call `ensure_pymol()` from your own app at startup:

    from install_pymol import ensure_pymol
    ensure_pymol()
"""

import sys
import subprocess
import importlib
import argparse


PACKAGE_NAME = "pymol-open-source"


def _try_import_pymol():
    """Try importing the PyMOL API. Returns (module_name, module) or (None, None)."""
    for mod_name in ("pymol2", "pymol"):
        try:
            mod = importlib.import_module(mod_name)
            return mod_name, mod
        except ImportError:
            continue
    return None, None


def _pip_install(package: str, upgrade: bool = False):
    """Run `<this interpreter> -m pip install <package>` as a subprocess."""
    cmd = [sys.executable, "-m", "pip", "install"]
    if upgrade:
        cmd.append("--upgrade")
    cmd.append(package)

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print("---- pip stdout ----")
        print(result.stdout)
        print("---- pip stderr ----")
        print(result.stderr)
        return False

    print(result.stdout)
    return True


def ensure_pymol(force: bool = False, check_only: bool = False) -> bool:
    """
    Ensure PyMOL is importable in this interpreter.

    Returns True if PyMOL is (or becomes) importable, False otherwise.
    """
    print(f"Target interpreter: {sys.executable}")
    print(f"Python version:     {sys.version.splitlines()[0]}")

    if not force:
        mod_name, mod = _try_import_pymol()
        if mod_name:
            path = getattr(mod, "__file__", "unknown location")
            print(f"✅ '{mod_name}' is already importable ({path})")
            return True
        else:
            print("❌ PyMOL not currently importable in this interpreter.")

    if check_only:
        print("Check-only mode: not installing.")
        return False

    print(f"Attempting: pip install {PACKAGE_NAME} into {sys.executable} ...")
    ok = _pip_install(PACKAGE_NAME, upgrade=force)

    if not ok:
        print(
            "\n⚠️  pip install failed.\n"
            "This is common on Windows, since PyMOL's open-source build has\n"
            "compiled C++/Qt dependencies that aren't always available as\n"
            "prebuilt pip wheels for every platform.\n\n"
            "Recommended fallback: use conda-forge instead, e.g.\n"
            "    conda create -n pymol-env python=3.11 -y\n"
            "    conda activate pymol-env\n"
            "    conda install -c conda-forge pymol-open-source -y\n"
        )
        return False

    # Re-check import after install
    mod_name, mod = _try_import_pymol()
    if mod_name:
        path = getattr(mod, "__file__", "unknown location")
        print(f"✅ Install succeeded. '{mod_name}' is now importable ({path})")
        return True

    print(
        "⚠️  pip reported success, but PyMOL still isn't importable.\n"
        "This can happen if pip installed into a different site-packages\n"
        "than the one this interpreter actually searches (unusual, but\n"
        "possible with unusual PYTHONPATH / sys.path setups). Try:\n"
        f"    {sys.executable} -c \"import pymol2\"\n"
        "and inspect the resulting error."
    )
    return False


def main():
    parser = argparse.ArgumentParser(description="Ensure PyMOL is installed in this interpreter.")
    parser.add_argument("--force", action="store_true", help="Reinstall/upgrade even if already present.")
    parser.add_argument("--check", action="store_true", help="Only check whether PyMOL is importable; don't install.")
    args = parser.parse_args()

    success = ensure_pymol(force=args.force, check_only=args.check)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()