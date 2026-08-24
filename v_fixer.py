from pathlib import Path
import re

# ===== EDIT THIS =====
ROOT = Path("ADD YOUR ROOT HERE")
# =====================

ses_pattern = re.compile(r"ses-(\d+)$")

for sub in ROOT.glob("sub-*"):
    if not sub.is_dir():
        continue

    for ses in sub.glob("ses-*"):
        name = ses.name

        # ses-00 -> ses-BL
        if name == "ses-00":
            new_name = "ses-BL"

        else:
            m = ses_pattern.match(name)
            if not m:
                continue

            num = int(m.group(1))
            new_name = f"ses-V{num:02d}"

        new_path = ses.parent / new_name

        print(f"{ses}  →  {new_path}")
        ses.rename(new_path)

print("Done renaming sessions.")