from __future__ import annotations
import argparse
import importlib
import sys
import glob
from pathlib import Path
from typing import Dict, List

# Map user-friendly tool names to package module paths
TOOLS: Dict[str, str] = {
    "gard": "EvolHarvest.gard",
    "busted": "EvolHarvest.busted",
    "fubar": "EvolHarvest.fubar",
    "fel": "EvolHarvest.fel",
    "meme": "EvolHarvest.meme",
    "absrel": "EvolHarvest.absrel",
}


def expand_inputs(pattern: str) -> List[str]:
    """Expand glob pattern; if no matches, return the pattern itself to let harvester handle it."""
    matches = glob.glob(pattern)
    return matches if matches else [pattern]


def run_tool(module_path: str, input_path: str, output_path: str, verbose: bool) -> None:
    """
    Load the module and call its standard `run(input_path, output_path, verbose=False)` function.
    Raises AttributeError if module doesn't expose run().
    """
    module = importlib.import_module(module_path)
    if not hasattr(module, "run"):
        raise AttributeError(
            f"Module {module_path} must define a `run(input_path, output_path, verbose=False)` function."
        )
    module.run(input_path, output_path, verbose=verbose)


def main(argv=None) -> int:
    """
    Command-line entry point for EvolHarvest.

    Usage examples:
      evolharvest --tool fel --input /path/to/jsons --output /path/to/outdir --verbose
      evolharvest fel --input /path/to/jsons --output /path/to/outdir --verbose
      evolharvest --list-tools
    """
    parser = argparse.ArgumentParser(prog="evolharvest", description="Run JSON harvesters for evolutionary analysis outputs.")
    # positional (optional) tool name
    parser.add_argument("tool_pos", nargs="?", default=None, help="Tool name to run (positional, alternative to --tool).")
    # named tool (backwards compatible)
    parser.add_argument("--tool", "-t", required=False, help="Tool name to run (or use --list-tools).")
    parser.add_argument("--input", "-i", required=False, help="Input file or glob pattern (e.g. 'results/*.json').")
    parser.add_argument("--output", "-o", required=False, help="Output file or directory. If multiple inputs, output must be a directory.")
    parser.add_argument("--list-tools", action="store_true", help="List available tool names.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    args = parser.parse_args(argv)

    # --list-tools handling
    if args.list_tools:
        print("Available tools:")
        for name in sorted(TOOLS.keys()):
            print(f"  {name} -> {TOOLS[name]}")
        return 0

    # determine tool: prefer positional, fall back to --tool
    tool = (args.tool_pos or args.tool)
    if not tool:
        parser.error("Tool name is required (provide as positional argument or with --tool).")
    tool = tool.lower()
    if tool not in TOOLS:
        parser.error(f"Unknown tool '{tool}'. Use --list-tools to see available names.")

    # validate input/output presence
    if not args.input:
        parser.error("--input is required")
    if not args.output:
        parser.error("--output is required")

    input_paths = expand_inputs(args.input)
    out_path = Path(args.output)

    # If multiple inputs, ensure output is (or will be) a directory
    if len(input_paths) > 1:
        # if the path exists and is not a dir -> error
        if out_path.exists() and not out_path.is_dir():
            parser.error("--output must be a directory when processing multiple input files")
        # if path doesn't exist, create it as a directory
        out_path.mkdir(parents=True, exist_ok=True)

        for inp in input_paths:
            inp_path = Path(inp)
            out_file = out_path / f"{inp_path.stem}_{tool}.csv"
            if args.verbose:
                print(f"[evolharvest] {tool}: {inp} -> {out_file}")
            # call harvester
            run_tool(TOOLS[tool], inp, str(out_file), verbose=args.verbose)

    else:
        # single input -> output can be a file or an existing directory
        if out_path.is_dir():
            inp_path = Path(input_paths[0])
            out_file = out_path / f"{inp_path.stem}_{tool}.csv"
        else:
            out_file = out_path
        if args.verbose:
            print(f"[evolharvest] {tool}: {input_paths[0]} -> {out_file}")
        run_tool(TOOLS[tool], input_paths[0], str(out_file), verbose=args.verbose)

    if args.verbose:
        print(f"[evolharvest] Completed {tool} on {len(input_paths)} input file(s).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
