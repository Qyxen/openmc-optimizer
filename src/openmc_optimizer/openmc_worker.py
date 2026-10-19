"""Run one prepared candidate model and calculate its PPF."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from openmc_optimizer.evaluator import evaluate_openmc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--statepoint", required=True)
    parser.add_argument("--assembly-cell-ids", required=True)
    parser.add_argument("--openmc-command", required=True)
    parser.add_argument("--ppf-tally-name", default="")
    parser.add_argument("--ppf-tally-score", default="fission-q-recoverable")
    args = parser.parse_args()

    command = json.loads(args.openmc_command)
    result = evaluate_openmc(
        args.directory.resolve(),
        args.statepoint,
        (
            [int(value) for value in args.assembly_cell_ids.split(",")]
            if args.assembly_cell_ids else None
        ),
        command,
        args.ppf_tally_name or None,
        args.ppf_tally_score,
    )
    (args.directory / "ppf_result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
