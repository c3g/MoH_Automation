#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def parse_mosdepth_metric(summary_path: Path) -> str:
    """
    Extract the 4th column of the last non-empty line from a mosdepth summary file.
    Args:
        summary_path (Path): Path to the mosdepth summary file.
    Returns:
        str: The extracted metric value.
    Raises:
        ValueError: If the file is empty or does not have the expected format.
    """
    with summary_path.open() as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        raise ValueError("Mosdepth summary file is empty")

    last_line = lines[-1]
    columns = last_line.split()

    if len(columns) < 4:
        raise ValueError(
            f"Unexpected mosdepth format: {last_line}"
        )

    return columns[3]


def construct_mosdepth_job(sample_name, job_start, job_stop, summary_path, metric):
    """
    Build the mosdepth job JSON object.
    Args:
        sample_name (str): Name of the sample.
        job_start (str): Job start time.
        job_stop (str): Job stop time.
        summary_path (Path): Path to the mosdepth summary file.
        metric (str): Extracted metric value.
    Returns:
        dict: The constructed mosdepth job JSON object.
    """
    summary_path = Path(summary_path)

    return {
        "job_name": f"mosdepth.{sample_name}",
        "job_start": job_start,
        "job_stop": job_stop,
        "job_status": "COMPLETED",
        "file": [
            {
                "location_uri": f"abacus://{summary_path}",
                "file_name": summary_path.name,
            }
        ],
        "metric": [
            {
                "metric_name": "dedup_coverage",
                "metric_value": metric,
            }
        ],
    }


def main():
    """
    Main function to add mosdepth job entries to tagged JSON.
    1. Parse command-line arguments.
    2. Extract metric from mosdepth summary file.
    3. Load tagged JSON file.
    4. Construct mosdepth job object.
    5. Traverse samples and readsets to append the job object.
    6. Write updated tagged JSON to a new file.
    """

    parser = argparse.ArgumentParser(
        description="Add mosdepth job entries to a tagged JSON file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Example usage:
      add_mosdepth_job.py \\    
        --sample-name SAMPLE_NAME \\
        --job-start "YYYY-MM-DD HH:MM:SS" \\
        --job-stop  "YYYY-MM-DD HH:MM:SS" \\
        --summary-txt /absolute/path/to/sample.mosdepth.summary.txt \\
        --tagged-json /absolute/path/to/tagged.json

    Notes:
      - The mosdepth metric is extracted from the last line, 4th column
        of the summary file.
      - All readsets whose name starts with the provided sample name
        will receive the new mosdepth job.
      - Output is written next to the input JSON with a '_mosdepth.json' suffix.
    """
    )

    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--job-start", required=True)
    parser.add_argument("--job-stop", required=True)
    parser.add_argument("--summary-txt", required=True)
    parser.add_argument("--tagged-json", required=True)

    args = parser.parse_args()

    summary_path = Path(args.summary_txt)
    tagged_json_path = Path(args.tagged_json)

    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)

    if not tagged_json_path.is_file():
        raise FileNotFoundError(tagged_json_path)

    # Parse metric
    metric_value = parse_mosdepth_metric(summary_path)

    # Load tagged JSON
    with tagged_json_path.open() as f:
        tagged = json.load(f)

    # Build job object template
    mosdepth_job = construct_mosdepth_job(
        sample_name=args.sample_name,
        job_start=args.job_start,
        job_stop=args.job_stop,
        summary_path=summary_path,
        metric=metric_value,
    )

    # Traverse samples / readsets
    for sample in tagged.get("sample", []):
        for readset in sample.get("readset", []):
            if readset.get("readset_name", "").startswith(args.sample_name):
                readset.setdefault("job", []).append(mosdepth_job)

    # Write output file
    output_path = tagged_json_path.with_name(
        tagged_json_path.stem + "_mosdepth.json"
    )
    with output_path.open("w") as f:
        json.dump(tagged, f, indent=4)


    print(f"Written: {output_path}")


if __name__ == "__main__":
    main()
