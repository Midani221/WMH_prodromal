from pathlib import Path
import pandas as pd
import re
import shutil
from collections import defaultdict

# ===== EDIT THESE =====
RAW_ROOTS = [
    Path(r"Roots of image data") # Change to the path containing the image data
]
CSV        = Path(r"CSV ROOT FOR COMPARISON") # Change the csv path to the file containing the subject IDs, and description names of the images
OUT        = Path(r"OUTPUT FILE") # Change to the output path  of interest
SUBJECT_ID = None  # set None for all subjects

PREFER_ROOT_INDEX = 0  # if tie, prefer RAW_ROOTS[0]
# ======================

OUT.mkdir(parents=True, exist_ok=True)

# ---------- CSV ----------
df = pd.read_csv(CSV).copy()

needed = ["Subject", "Visit", "Description", "Acq Date", "T1", "FLAIR"]
missing = [c for c in needed if c not in df.columns]
assert not missing, f"Missing columns: {missing}"

df["Subject"] = df["Subject"].astype(str).str.strip()
df["Description"] = df["Description"].astype(str).str.strip()

# parse date; if your CSV is m/d/Y, you can speed+stabilize with format="%m/%d/%Y"
df["Acq Date"] = pd.to_datetime(df["Acq Date"], errors="coerce")
df["acq_date_key"] = df["Acq Date"].dt.strftime("%Y-%m-%d")

if SUBJECT_ID is not None:
    df = df[df["Subject"] == str(SUBJECT_ID)].copy()
    assert len(df) > 0, f"No rows found for Subject={SUBJECT_ID}"

# ---------- Helpers ----------
dt_folder_re = re.compile(r"(?P<date>\d{4}-\d{2}-\d{2})_(?P<time>\d{2}_\d{2}_\d{2})(?:\.\d+)?$")

def folder_date_key(p: Path):
    m = dt_folder_re.match(p.name)
    return m.group("date") if m else None

def ses_from_visit(v):
    if pd.isna(v):
        return "ses-00"
    m = re.search(r"(\d+)", str(v))
    return f"ses-{int(m.group(1)):02d}" if m else "ses-00"

def truthy(x):
    if pd.isna(x): return False
    if isinstance(x, (int, float)): return x != 0
    return str(x).strip().lower() not in ("0","no","n","false","")

def norm_sub(subject):
    s = str(subject).strip()
    return s if s.startswith("sub-") else f"sub-{s}"

def count_files(folder: Path):
    return sum(1 for _ in folder.rglob("*") if _.is_file())

# Optional: enable if you want to tolerate 3D SAG vs 3D_SAG etc.
# def norm_desc(s: str) -> str:
#     return re.sub(r"[^a-z0-9]+", "", str(s).lower())
# USE_DESC_NORMALIZATION = True

USE_DESC_NORMALIZATION = False
def desc_key(s: str) -> str:
    if not USE_DESC_NORMALIZATION:
        return str(s).strip()
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())

# ---------- Build disk index ONCE ----------
# index[(subject, desc_key, date_key)] -> list of (root_index, datetime_folder_path)
disk_index = defaultdict(list)

def build_disk_index():
    # If you're running per subject, this speeds further by only scanning that subject folder:
    subjects_to_scan = set(df["Subject"].unique())

    for i, root in enumerate(RAW_ROOTS):
        for subject in subjects_to_scan:
            subroot = root / subject
            if not subroot.exists():
                continue

            # Find datetime folders anywhere under subject
            for p in subroot.rglob("*"):
                if not p.is_dir():
                    continue
                dkey = folder_date_key(p)
                if not dkey:
                    continue
                # Expect: .../<subject>/<description>/<datetimefolder>
                desc = p.parent.name
                disk_index[(subject, desc_key(desc), dkey)].append((i, p))

    # deterministic ordering
    for k in list(disk_index.keys()):
        disk_index[k].sort(key=lambda x: (x[1].name, x[0]))

build_disk_index()

def find_hits(subject: str, description: str, acq_date_key: str):
    return disk_index.get((subject, desc_key(description), acq_date_key), [])

def choose_best_hit(hits):
    scored = []
    for i, p in hits:
        scored.append((count_files(p), 1 if i == PREFER_ROOT_INDEX else 0, i, p))
    scored.sort(reverse=True)
    best = scored[0]
    return best[3], best[2]  # (path, root_index)

# ---------- MOVE ----------
issues = []
moved_sources = set()
staged = []

for r in df.itertuples(index=False):
    subject = str(r.Subject).strip()
    desc = str(r.Description).strip()
    ses = ses_from_visit(r.Visit)
    date_key = r.acq_date_key

    if pd.isna(date_key):
        issues.append(("bad_csv_date", subject, desc))
        continue

    hits = find_hits(subject, desc, date_key)
    if not hits:
        issues.append(("not_found_in_any_root", subject, ses, desc, date_key))
        continue

    if len(hits) > 1:
        issues.append(("multiple_hits", subject, ses, desc, date_key, [(i, str(p)) for i, p in hits]))

    src, src_root_idx = choose_best_hit(hits)

    if str(src) in moved_sources:
        issues.append(("already_moved_skip", subject, ses, desc, date_key, str(src)))
        continue

    mods = []
    if truthy(r.T1): mods.append("T1")
    if truthy(r.FLAIR): mods.append("FLAIR")
    if not mods: mods = ["UNKNOWN"]

    if len(mods) > 1:
        issues.append(("row_marks_multiple_modalities", subject, ses, desc, date_key, mods))
        mods = ["RAW"]

    mod = mods[0]
    dst_root = OUT / norm_sub(subject) / ses / "DICOM" / mod / desc
    dst_root.mkdir(parents=True, exist_ok=True)

    dst = dst_root / src.name
    if dst.exists():
        issues.append(("dest_exists_skip", str(dst)))
        continue

    shutil.move(str(src), str(dst))
    moved_sources.add(str(src))
    staged.append((subject, ses, mod, desc, date_key, src_root_idx, str(dst)))

print("Moved folders:", len(staged))
print("Issues flagged:", len(issues))
print("Output root:", OUT)

# Optional: peek at first issues
from pprint import pprint
pprint(issues[:30])