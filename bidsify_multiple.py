import argparse
import csv
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Tuple, Optional, List, Set

def require_cmd(cmd: str) -> None:
    if shutil.which(cmd) is None:
        raise RuntimeError(f"Command not found on PATH: {cmd}. Activate your dcm2bids_env.")

def run_cmd(cmd, cwd: Optional[Path] = None) -> None:
    cmd = [str(x) for x in cmd]
    print("\n>>>", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd) if cwd else None)
    if res.returncode != 0:
        print("!!! FAILED:", res.returncode)
        print("--- STDOUT (first 8000 chars) ---")
        print(res.stdout[:8000])
        print("--- STDERR (first 8000 chars) ---")
        print(res.stderr[:8000])
        raise subprocess.CalledProcessError(res.returncode, cmd, output=res.stdout, stderr=res.stderr)

def is_yes(x: str) -> bool:
    return str(x).strip().lower() in {"yes", "y", "true", "1"}

def visit_to_session_label(visit: str) -> str:
    """
    Convert CSV Visit into dcm2bids session label (without 'ses-').
    Examples:
      "00", "0", "BL", "baseline" -> "BL"
      "6", "06", "V06", "ses-V06" -> "V06"
    """
    v = str(visit).strip()
    v_low = v.lower()

    if v_low in {"bl", "baseline", "ses-bl"}:
        return "BL"

    # if already ses-xxx
    if v_low.startswith("ses-"):
        v = v[4:]

    # if already Vxx
    if re.fullmatch(r"v\d{2}", v.lower()):
        return v.upper()

    # numeric -> Vxx
    m = re.search(r"(\d+)", v)
    if m:
        num = int(m.group(1))
        if num == 0:
            return "BL"
        return f"V{num:02d}"

    return "BL"

def write_config(config_path: Path, t1_descs: List[str], flair_descs: List[str]) -> None:
    """
    Robust config: multiple possible SeriesDescriptions per suffix.
    dcm2bids will match whichever exists in that session input.
    """
    descriptions = []
    for d in t1_descs:
        descriptions.append({"datatype": "anat", "suffix": "T1w", "criteria": {"SeriesDescription": d}})
    for d in flair_descs:
        descriptions.append({"datatype": "anat", "suffix": "FLAIR", "criteria": {"SeriesDescription": d}})

    cfg = {"descriptions": descriptions}
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    print("Config saved:", config_path)

def main():
    ap = argparse.ArgumentParser(
        description="BIDSify PPMI subset from CSV, but run dcm2bids PER SESSION using staged layout sub-*/ses-*/DICOM to avoid mixing."
    )
    ap.add_argument("--project_root", required=True, help="Working folder (configs/ and bids/ will be created here).")
    ap.add_argument("--ppmi_root", required=True, help="Staged root containing sub-*/ses-*/DICOM (NOT the raw PPMI download root).")
    ap.add_argument("--csv", required=True, help="Subset CSV with headers: Subject, Visit, Description, T1, FLAIR.")
    ap.add_argument("--dry_run", default="0", help="1 = print commands only, no conversion.")
    args = ap.parse_args()

    project_root = Path(args.project_root).resolve()
    ppmi_root = Path(args.ppmi_root).resolve()
    csv_path = Path(args.csv).resolve()
    dry_run = str(args.dry_run).strip().lower() in {"1", "yes", "true", "y"}

    for c in ["dcm2bids", "dcm2bids_scaffold"]:
        require_cmd(c)

    bids_out = project_root / "bids"
    bids_out.mkdir(parents=True, exist_ok=True)

    # Scaffold once
    scaffold_cmd = ["dcm2bids_scaffold", "-o", str(bids_out)]
    if dry_run:
        print("[DRY RUN] Would scaffold:", " ".join(scaffold_cmd))
    else:
        try:
            run_cmd(scaffold_cmd, cwd=project_root)
        except subprocess.CalledProcessError:
            if (bids_out / "dataset_description.json").exists():
                print("Scaffold failed but dataset_description.json exists; continuing.")
            else:
                raise

    required_cols = {"Subject", "Visit", "Description", "T1", "FLAIR"}

    # Group rows per (Subject, Visit) but store *sets* of possible SeriesDescriptions
    groups: Dict[Tuple[str, str], Dict[str, Set[str]]] = {}

    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError("CSV has no headers.")
        missing = required_cols - set(reader.fieldnames)
        if missing:
            raise ValueError(f"CSV missing required columns: {sorted(missing)}")

        for row in reader:
            sub = str(row["Subject"]).strip()
            visit = str(row["Visit"]).strip()
            desc = str(row["Description"]).strip()

            if not sub or not visit or not desc:
                continue

            key = (sub, visit)
            g = groups.setdefault(key, {"t1_descs": set(), "flair_descs": set()})

            if is_yes(row["T1"]):
                g["t1_descs"].add(desc)
            if is_yes(row["FLAIR"]):
                g["flair_descs"].add(desc)

    print(f"\nGroups found (Subject+Visit): {len(groups)}")

    for (sub, visit), g in sorted(groups.items()):
        t1_descs = sorted(g["t1_descs"])
        flair_descs = sorted(g["flair_descs"])
        ses_label = visit_to_session_label(visit)  # e.g., BL or V06

        if not t1_descs or not flair_descs:
            print(f"Skipping Subject={sub}, Visit={visit}: missing T1 or FLAIR mapping in CSV.")
            continue

        # staged layout root:
        # ppmi_root/sub-<sub>/ses-<ses_label>/DICOM
        sub_dir = ppmi_root / f"sub-{sub}"
        if not sub_dir.exists():
            # fallback if your sub folders are not prefixed
            alt = ppmi_root / sub
            if alt.exists():
                sub_dir = alt
            else:
                print(f"Skipping {sub}: subject folder not found under {ppmi_root}")
                continue

        ses_dir = sub_dir / f"ses-{ses_label}"
        dicom_dir = ses_dir / "DICOM"
        if not dicom_dir.exists():
            print(f"Skipping sub-{sub} ses-{ses_label}: DICOM folder not found: {dicom_dir}")
            continue

        print("\n" + "=" * 90)
        print(f"Processing Subject={sub}  Visit='{visit}'  SessionLabel={ses_label}")
        print("DICOM input root:", dicom_dir)
        print("T1 SeriesDescriptions   :", t1_descs)
        print("FLAIR SeriesDescriptions:", flair_descs)

        config_path = project_root / "configs" / f"sub-{sub}_ses-{ses_label}_config.json"
        write_config(config_path, t1_descs, flair_descs)

        cmd = [
            "dcm2bids",
            "-d", str(dicom_dir),      # IMPORTANT: session-isolated input
            "-p", str(sub),
            "-s", str(ses_label),      # IMPORTANT: pass BL or V06 (NOT 06)
            "-c", str(config_path),
            "-o", str(bids_out),
        ]

        if dry_run:
            print("[DRY RUN] Would run:", " ".join(cmd))
        else:
            run_cmd(cmd, cwd=project_root)
            print(f"✅ Done sub-{sub} ses-{ses_label}")

    print("\nAll done. Output:", bids_out)

if __name__ == "__main__":
    main()