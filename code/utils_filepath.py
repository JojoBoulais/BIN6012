from pathlib import Path
import os


CURRENT_DIR = Path(__file__).parent

KINDS = ["grn", "data", "model", "tmp", "results"]

def get_filepath(kind, name=""):

    if kind not in KINDS:
        print(f"kind must be one of the following: " + ",".join(KINDS))
        return ""

    return os.path.join(CURRENT_DIR, "..", kind, name)


