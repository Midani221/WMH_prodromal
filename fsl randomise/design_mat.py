#!/usr/bin/env python3

import csv
from pathlib import Path
from statistics import mean


# ============================================================
# INPUTS
# ============================================================

INPUT_TABLE = Path("Design_matrix.csv")
INPUT_IMAGE_LIST = Path("image_list.txt")
OUTPUT_ROOT = Path("designs")

SUBJECT_COLUMN = "SubjID"
AGE_COLUMN = "baseline_age"
SEX_COLUMN = "Sx"

# Predictor -> corresponding baseline variable
MODELS = {
    "phenoconversion": {
        "predictor": "phenoconversion_binary",
        "baseline": None,
        "binary_predictor": True,
    },

    "FPC1_z": {
        "predictor": "FPC1_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC2_z": {
        "predictor": "FPC2_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC_scopa1_z": {
        "predictor": "FPC_scopa1_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC_scopa2_z": {
        "predictor": "FPC_scopa2_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC_moca1_z": {
        "predictor": "FPC_moca1_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC_moca2_z": {
        "predictor": "FPC_moca2_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC_EF1_z": {
        "predictor": "FPC_EF1_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC_EF2_z": {
        "predictor": "FPC_EF2_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC_M1_z": {
        "predictor": "FPC_M1_z",
        "baseline": None,
        "binary_predictor": False,
    },

    "FPC_M2_z": {
        "predictor": "FPC_M2_z",
        "baseline": None,
        "binary_predictor": False,
    },
}


MISSING_VALUES = {
    "",
    "na",
    "n/a",
    "nan",
    "null",
    ".",
}


def parse_number(value):
    """Convert a value to float, returning None for missing values."""
    value = str(value).strip()

    if value.lower() in MISSING_VALUES:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def read_table(path):
    """Read comma-separated or tab-separated files."""
    with path.open("r", newline="", encoding="utf-8-sig") as file:
        sample = file.read(4096)
        file.seek(0)

        try:
            dialect = csv.Sniffer().sniff(
                sample,
                delimiters=",\t;"
            )
        except csv.Error:
            dialect = csv.excel

        reader = csv.DictReader(file, dialect=dialect)

        if reader.fieldnames is None:
            raise ValueError("Input table has no header.")

        return list(reader.fieldnames), list(reader)


# ============================================================
# READ TABLE AND IMAGE PATHS
# ============================================================

columns, rows = read_table(INPUT_TABLE)

with INPUT_IMAGE_LIST.open(
    "r",
    encoding="utf-8-sig"
) as image_file:

    image_paths = [
        line.strip().rstrip("\r")
        for line in image_file
        if line.strip()
    ]


if len(rows) != len(image_paths):
    raise ValueError(
        "CSV rows and image paths do not match.\n"
        f"CSV rows: {len(rows)}\n"
        f"Image paths: {len(image_paths)}"
    )


if SUBJECT_COLUMN not in columns:
    raise ValueError(
        f"Subject column '{SUBJECT_COLUMN}' was not found."
    )


# Attach each image to its corresponding row
for row, image_path in zip(rows, image_paths):
    row["_image_path"] = image_path


OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)


# ============================================================
# GENERATE EACH MODEL
# ============================================================

generated_models = []
skipped_models = []

for model_name, settings in MODELS.items():

    predictor = settings["predictor"]
    baseline = settings["baseline"]
    binary_predictor = settings["binary_predictor"]

    model_variables = [
        predictor,
        AGE_COLUMN,
        SEX_COLUMN,
    ]

    if baseline is not None:
        model_variables.insert(1, baseline)

    required_columns = set(model_variables)
    missing_columns = required_columns.difference(columns)

    if missing_columns:
        print(
            f"SKIPPING {model_name}: "
            f"missing columns {sorted(missing_columns)}"
        )
        skipped_models.append(model_name)
        continue

    complete_rows = []

    for row in rows:

        parsed = {}
        complete = True

        for variable in model_variables:
            value = parse_number(row.get(variable, ""))

            if value is None:
                complete = False
                break

            parsed[variable] = value

        if not complete:
            continue

        # Sex must be 0/1
        if parsed[SEX_COLUMN] not in (0, 1):
            raise ValueError(
                f"{SEX_COLUMN} must be 0 or 1. "
                f"Subject: {row[SUBJECT_COLUMN]}, "
                f"value: {parsed[SEX_COLUMN]}"
            )

        # Phenoconversion must be 0/1
        if binary_predictor:
            if parsed[predictor] not in (0, 1):
                raise ValueError(
                    f"{predictor} must be 0 or 1. "
                    f"Subject: {row[SUBJECT_COLUMN]}, "
                    f"value: {parsed[predictor]}"
                )

        row_copy = dict(row)
        row_copy["_parsed"] = parsed
        complete_rows.append(row_copy)

    if not complete_rows:
        print(f"SKIPPING {model_name}: no complete observations")
        skipped_models.append(model_name)
        continue

    predictor_values = [
        row["_parsed"][predictor]
        for row in complete_rows
    ]

    if len(set(predictor_values)) < 2:
        print(
            f"SKIPPING {model_name}: "
            "predictor has no variation"
        )
        skipped_models.append(model_name)
        continue

    # Mean-centre age and corresponding baseline measurement
    variables_to_center = [AGE_COLUMN]

    if baseline is not None:
        variables_to_center.append(baseline)

    centering_means = {
        variable: mean(
            row["_parsed"][variable]
            for row in complete_rows
        )
        for variable in variables_to_center
    }

    # Design order:
    # Intercept, predictor, baseline if present, age, sex
    design_columns = ["Intercept", predictor]

    if baseline is not None:
        design_columns.append(f"{baseline}_centered")

    design_columns.extend([
        f"{AGE_COLUMN}_centered",
        SEX_COLUMN,
    ])

    model_directory = OUTPUT_ROOT / model_name
    model_directory.mkdir(parents=True, exist_ok=True)

    design_path = model_directory / "design.txt"
    contrast_path = model_directory / "contrasts.txt"
    subject_path = model_directory / "subjects.csv"
    image_list_path = model_directory / "image_list.txt"
    column_path = model_directory / "design_columns.txt"

    # --------------------------------------------------------
    # Design matrix
    # --------------------------------------------------------

    with design_path.open("w", encoding="utf-8") as file:

        for row in complete_rows:

            parsed = row["_parsed"]

            design_values = [
                1.0,
                parsed[predictor],
            ]

            if baseline is not None:
                design_values.append(
                    parsed[baseline]
                    - centering_means[baseline]
                )

            design_values.extend([
                parsed[AGE_COLUMN]
                - centering_means[AGE_COLUMN],

                parsed[SEX_COLUMN],
            ])

            file.write(
                " ".join(
                    f"{value:.10g}"
                    for value in design_values
                )
                + "\n"
            )

    # --------------------------------------------------------
    # Positive and negative predictor contrasts
    # --------------------------------------------------------

    number_of_evs = len(design_columns)

    positive_contrast = [0] * number_of_evs
    negative_contrast = [0] * number_of_evs

    # Predictor is always EV2
    positive_contrast[1] = 1
    negative_contrast[1] = -1

    with contrast_path.open("w", encoding="utf-8") as file:
        file.write(
            " ".join(map(str, positive_contrast)) + "\n"
        )
        file.write(
            " ".join(map(str, negative_contrast)) + "\n"
        )

    # --------------------------------------------------------
    # Save column definitions
    # --------------------------------------------------------

    with column_path.open("w", encoding="utf-8") as file:

        for number, column in enumerate(
            design_columns,
            start=1
        ):
            file.write(f"EV{number}: {column}\n")

        file.write("\n")
        file.write(
            "Contrast 1: positive association/"
            "converter greater than non-converter\n"
        )
        file.write(
            "Contrast 2: negative association/"
            "non-converter greater than converter\n"
        )

        file.write("\nCentring means:\n")

        for variable, value in centering_means.items():
            file.write(f"{variable}: {value:.6f}\n")

    # --------------------------------------------------------
    # Model-specific subject CSV
    # --------------------------------------------------------

    with subject_path.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        output_columns = columns + ["image_path"]

        writer = csv.DictWriter(
            file,
            fieldnames=output_columns
        )

        writer.writeheader()

        for row in complete_rows:
            output_row = {
                column: row.get(column, "")
                for column in columns
            }

            output_row["image_path"] = row["_image_path"]
            writer.writerow(output_row)

    # --------------------------------------------------------
    # Model-specific image list
    # --------------------------------------------------------

    with image_list_path.open(
        "w",
        encoding="utf-8",
        newline="\n"
    ) as file:

        for row in complete_rows:
            file.write(row["_image_path"] + "\n")

    generated_models.append(model_name)

    print(
        f"GENERATED {model_name}: "
        f"{len(complete_rows)} participants, "
        f"{len(design_columns)} EVs"
    )


# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("Generated models:")
for model in generated_models:
    print(f"  {model}")

if skipped_models:
    print()
    print("Skipped models:")
    for model in skipped_models:
        print(f"  {model}")