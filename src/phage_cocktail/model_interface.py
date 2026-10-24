"""Command-line entry point for the ModelInputBundle / PairPrediction v1 exchange."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .artifacts import atomic_write_json, read_json_object
from .input_bundle import build_pair_prediction, load_input_bundle, validate_pair_prediction


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phage-model-interface",
        description="Validate exchange files and emit explicit unresolved pair predictions.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    bundle = commands.add_parser(
        "validate-bundle", help="verify frozen snapshots and every local payload"
    )
    bundle.add_argument("bundle", type=Path)
    bundle.add_argument(
        "--require-sequence", choices=("host", "phage"), action="append", default=[]
    )

    pair = commands.add_parser(
        "validate-pair-prediction", help="check an exact source-endpoint result"
    )
    pair.add_argument("result", type=Path)
    pair.add_argument("--bundle", type=Path, required=True)

    abstain = commands.add_parser(
        "create-unresolved-pair-prediction", help="record no inference against a frozen bundle"
    )
    abstain.add_argument("bundle", type=Path)
    abstain.add_argument("output", type=Path)
    abstain.add_argument("--prediction-id", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        bundle = load_input_bundle(
            args.bundle, required_sequences=getattr(args, "require_sequence", ())
        )
        if args.command == "validate-bundle":
            print(f"valid ModelInputBundle: {bundle['bundle_id']} sha256={bundle['content_sha256']}")
        elif args.command == "validate-pair-prediction":
            result = validate_pair_prediction(read_json_object(args.result), bundle)
            print(f"valid PairPrediction: {result['prediction_id']}; {result['estimate']['status']}")
        else:
            if args.output.exists() or args.output.is_symlink():
                raise ValueError("output already exists; refusing to overwrite")
            if args.output.resolve().is_relative_to(args.bundle.resolve().parent):
                raise ValueError("write prediction outside the frozen bundle directory")
            result = build_pair_prediction(bundle, prediction_id=args.prediction_id)
            atomic_write_json(result, args.output)
            print(f"wrote unresolved PairPrediction: {args.output}; no inference was run")
    except (ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
