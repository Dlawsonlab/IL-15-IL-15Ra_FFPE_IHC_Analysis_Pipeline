"""
Project paths for the Python scripts.

ROOT is the project folder, found from this file's location. Site-specific paths come from
config/site.env, the same file the hpc/ wrappers source. Decisions the pipeline reads as
input live in data/decisions/.

Scripts in the role folders under scripts/ put scripts/ on sys.path to import this module.
"""
import csv
import os
import shlex

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MACROS = os.path.join(ROOT, "macros")
DECISIONS = os.path.join(ROOT, "data", "decisions")
FIELD_EXCLUSIONS = os.path.join(DECISIONS, "field_exclusions")


def _read_site_env():
    site = {}
    with open(os.path.join(ROOT, "config", "site.env"), encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            words = shlex.split(value, comments=True)
            site[key.strip()] = words[0] if words else ""
    return site


SITE = _read_site_env()


def to_local(path):
    """Map a cluster path from imagelist.txt onto the workstation's mount of the same share."""
    return path.replace(SITE["HPC_SHARE"], SITE["LOCAL_SHARE"], 1)


def field_exclusions(sub):
    """Field-exclusion file for a per_sample_out/ subfolder. File names drop any _EXCLUDED suffix."""
    return os.path.join(FIELD_EXCLUSIONS, f"{sub.removesuffix('_EXCLUDED')}_fov_exclude.csv")


def field_exclusion_rows(sub):
    """A section's excluded fields as rows (fov_id, reason); empty if it has none."""
    path = field_exclusions(sub)
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))
