from pathlib import Path
import pandas as pd

# ============================================================
# SETTINGS
# ============================================================

CSV_FILE = Path("") #Subjects file containing id and ses

# Folders containing the subject directories.
# Add as many search locations as needed.
SEARCH_ROOTS = [


]

# Change these to match your CSV column names
SUBJECT_COLUMN = "Sub"
SESSION_COLUMN = "Ses"

TARGET_FILENAME = "results2mni_combined_nonlin.nii.gz"

OUTPUT_CSV = Path("subjects_with_image_paths.csv")
MATCHED_CSV = Path("matched_subjects_only.csv")
IMAGE_LIST = Path("image_list.txt")


# ============================================================
# FUNCTIONS
# ============================================================

def add_prefix(value, prefix):
    """Add sub- or ses- if the prefix is not already present."""
    value = str(value).strip()

    if not value:
        return ""

    if value.lower().startswith(prefix.lower()):
        return value

    return f"{prefix}{value}"


# ============================================================
# READ CSV
# ============================================================

# dtype=str prevents subject numbers from being converted to numbers
df = pd.read_csv(
    CSV_FILE,
    dtype=str,
    keep_default_na=False
)

required_columns = {SUBJECT_COLUMN, SESSION_COLUMN}

missing_columns = required_columns.difference(df.columns)

if missing_columns:
    raise ValueError(
        f"Missing CSV columns: {sorted(missing_columns)}\n"
        f"Available columns: {list(df.columns)}"
    )

df["subject_folder"] = df[SUBJECT_COLUMN].apply(
    lambda x: add_prefix(x, "sub-")
)

df["session_folder"] = df[SESSION_COLUMN].apply(
    lambda x: add_prefix(x, "ses-")
)


# ============================================================
# INDEX ALL TARGET FILES
# ============================================================

file_index = {}

for root in SEARCH_ROOTS:

    if not root.exists():
        print(f"WARNING: Search folder does not exist: {root}")
        continue

    print(f"Searching: {root}")

    for path in root.rglob(TARGET_FILENAME):

        # Expected structure:
        # sub-xxxx/ses-xxxx/output/results2mni_nonlin.nii.gz

        if path.parent.name != "output":
            continue

        session = path.parent.parent.name
        subject = path.parent.parent.parent.name

        if not subject.startswith("sub-"):
            continue

        if not session.startswith("ses-"):
            continue

        key = (subject.lower(), session.lower())

        file_index.setdefault(key, []).append(path.resolve())


# ============================================================
# MATCH CSV ROWS TO FILES
# ============================================================

image_paths = []
statuses = []
number_of_matches = []

for _, row in df.iterrows():

    subject = row["subject_folder"]
    session = row["session_folder"]

    key = (subject.lower(), session.lower())

    # Remove duplicated representations of the same path
    matches = sorted(set(file_index.get(key, [])))

    number_of_matches.append(len(matches))

    if len(matches) == 0:
        image_paths.append("")
        statuses.append("missing")

    elif len(matches) == 1:
        image_paths.append(str(matches[0]))
        statuses.append("found")

    else:
        # Do not silently choose between duplicate files
        image_paths.append(" | ".join(str(x) for x in matches))
        statuses.append("multiple_matches")


df["image_path"] = image_paths
df["match_status"] = statuses
df["number_of_matches"] = number_of_matches


# ============================================================
# SAVE RESULTS
# ============================================================

# Complete report, including missing and duplicate matches
df.to_csv(OUTPUT_CSV, index=False)

# Only unambiguous matches, preserving original CSV order
matched_df = df.loc[df["match_status"] == "found"].copy()

matched_df.to_csv(MATCHED_CSV, index=False)

# Plain-text image list for fslmerge
matched_df["image_path"].to_csv(
    IMAGE_LIST,
    index=False,
    header=False
)


# ============================================================
# SUMMARY
# ============================================================

print("\nSearch complete")
print("----------------")
print(f"CSV rows:          {len(df)}")
print(f"Files found:       {(df['match_status'] == 'found').sum()}")
print(f"Missing files:     {(df['match_status'] == 'missing').sum()}")
print(f"Multiple matches:  {(df['match_status'] == 'multiple_matches').sum()}")

print(f"\nFull report:       {OUTPUT_CSV.resolve()}")
print(f"Matched rows:      {MATCHED_CSV.resolve()}")
print(f"FSL image list:    {IMAGE_LIST.resolve()}")