from pathlib import Path
from warnings import warn
from hashlib import file_digest, sha256

import nbformat
from nbconvert.exporters import ScriptExporter

_HASH_TYPE = "sha256"


def _hash_of_file(file: Path) -> str:
    with file.open("rb") as f:
        return file_digest(f, _HASH_TYPE).hexdigest()


def _read_hash(file: Path) -> str:
    if not file.is_file():
        return ""
    else:
        return file.read_text().strip()


def generate_py_from_ipynb(app) -> None:
    _base_dir = Path(app.srcdir, "contents")
    notebooks_dir = _base_dir / "notebooks"
    scripts_dir = _base_dir / "scripts"

    if not notebooks_dir.is_dir():
        warn(f"[MyST-NB Export] Skipped: notebooks directory not found: {notebooks_dir}")
        return

    scripts_dir.mkdir(exist_ok=True)

    exporter = ScriptExporter()

    for ipynb in notebooks_dir.glob("*.ipynb"):
        # > Determine path to output script
        stem = ipynb.stem
        script_path = scripts_dir / f"{stem}.py"

        # > Calculate hash of notebook
        ipynb_hash = _hash_of_file(ipynb)
        hash_file = script_path.with_suffix(f".{_HASH_TYPE}")

        # > Read hash (if exists) and check if it changed
        read_hash = _read_hash(hash_file)
        # > If both hashes are equal skip conversion
        if read_hash and read_hash == ipynb_hash:
            continue
        else:
            # > Delete script and hash file.
            script_path.unlink(missing_ok=True)
            hash_file.unlink(missing_ok=True)

            # > Convert Jupyter notebook to a script
            with ipynb.open() as f_ipynb:
                nb_node = nbformat.read(f_ipynb, as_version=4)
            # << END OF WITH

            (body, _) = exporter.from_notebook_node(nb_node)

            # > Dump script and hash
            script_path.write_text(body)
            hash_file.write_text(ipynb_hash)

            print(f"[MyST-NB Export] Exported: {script_path} -> scripts/{stem}.py")
