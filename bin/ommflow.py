#!/usr/bin/env python3
"""Command-line entry point for the OpenMM protein MD workflow."""

from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from ommflow.lib.config import parse_arguments, write_default_configuration
from ommflow.lib.simulation import run_workflow
from ommflow.lib.system_builder import describe_input_components


def main() -> None:
    """Parse command-line options and run the requested workflow."""
    args = parse_arguments()
    if args.write_default_config:
        write_default_configuration(args.write_default_config)
        print(f"Default configuration written to {args.write_default_config}")
        return
    if args.list_components:
        print(describe_input_components(args))
        return
    run_workflow(args)


if __name__ == "__main__":
    main()
