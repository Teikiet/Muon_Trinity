#!/usr/bin/env python3

"""Run clang-tidy with the CORSIKA policy and a CMake compilation database.

By default, check C++ translation units changed relative to ``main``. Use
``--all`` in continuous integration to check every configured translation unit.
The build directory must contain ``compile_commands.json``.
"""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


SOURCE_SUFFIXES = {".cc", ".cpp", ".cxx", ".c++"}
HEADER_SUFFIXES = {".h", ".hh", ".hpp", ".hxx", ".inl"}
EXCLUDED_DIRECTORIES = {".git", "build", "externals", "install"}


def parse_arguments(repo_root: Path) -> argparse.Namespace:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
      "--all",
      action="store_true",
      help="check all translation units from the compilation database",
  )
  parser.add_argument(
      "--base-ref",
      default="main",
      help="Git revision used to select changed files (default: main)",
  )
  parser.add_argument(
      "--build-dir",
      type=Path,
      default=repo_root / "build",
      help="CMake build directory containing compile_commands.json",
  )
  parser.add_argument(
      "--fix",
      action="store_true",
      help="apply clang-tidy fixes; review the resulting changes before committing",
  )
  parser.add_argument(
      "--header-filter",
      help="regular expression selecting project headers to diagnose",
  )
  return parser.parse_args()


def is_project_source(path: Path, repo_root: Path) -> bool:
  try:
    relative_path = path.resolve().relative_to(repo_root)
  except ValueError:
    return False
  return not any(part in EXCLUDED_DIRECTORIES for part in relative_path.parts)


def project_source_root(path: Path, repo_root: Path) -> Path | None:
  """Return the checkout root that contains a project source path.

  Compilation databases can be generated through a different mount point than
  the one used to run this script. Keep that database path for clang-tidy, but
  recognize it when its path below the repository directory names a source in
  this checkout.
  """
  if is_project_source(path, repo_root):
    return repo_root

  for index, part in enumerate(path.parts):
    if part != repo_root.name:
      continue
    relative_path = Path(*path.parts[index + 1 :])
    candidate = repo_root / relative_path
    if candidate.is_file() and is_project_source(candidate, repo_root):
      return Path(*path.parts[: index + 1])
  return None


def project_relative_path(path: Path, repo_root: Path) -> Path | None:
  source_root = project_source_root(path, repo_root)
  if source_root is None:
    return None
  return path.relative_to(source_root)


def changed_files(repo_root: Path, base_ref: str) -> list[Path]:
  command = ["git", "diff", "--name-only", "--diff-filter=ACMR", base_ref]
  result = subprocess.run(
      command,
      cwd=repo_root,
      check=True,
      capture_output=True,
      text=True,
  )
  files = [repo_root / line for line in result.stdout.splitlines()]

  untracked = subprocess.run(
      ["git", "ls-files", "--exclude-standard", "--others"],
      cwd=repo_root,
      check=True,
      capture_output=True,
      text=True,
  )
  files.extend(repo_root / line for line in untracked.stdout.splitlines())
  return files


def compilation_files(database_path: Path, repo_root: Path) -> list[Path]:
  try:
    entries = json.loads(database_path.read_text(encoding="utf-8"))
  except json.JSONDecodeError as error:
    raise RuntimeError(f"invalid compilation database: {error}") from error

  files = set()
  for entry in entries:
    directory = Path(entry["directory"])
    source_file = Path(entry["file"])
    if not source_file.is_absolute():
      source_file = directory / source_file
    source_file = source_file.resolve()
    if (
        source_file.suffix in SOURCE_SUFFIXES
        and project_source_root(source_file, repo_root) is not None
    ):
      files.add(source_file)
  return sorted(files)


def default_header_filter(project_roots: set[Path]) -> str:
  source_directories = "|".join(
      re.escape(str(project_root / directory))
      for project_root in sorted(project_roots)
      for directory in ("applications", "corsika", "examples", "src", "tests")
  )
  return f"^({source_directories})/"


def main() -> int:
  repo_root = Path(__file__).resolve().parent
  arguments = parse_arguments(repo_root)
  build_dir = arguments.build_dir.resolve()
  database_path = build_dir / "compile_commands.json"
  if not database_path.is_file():
    print(
        f"Compilation database not found: {database_path}\n"
        "Configure CMake with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON or pass --build-dir.",
        file=sys.stderr,
    )
    return 2

  clang_tidy = os.environ.get("CLANG_TIDY", "clang-tidy")
  if shutil.which(clang_tidy) is None:
    print(f"clang-tidy executable not found: {clang_tidy}", file=sys.stderr)
    return 2

  try:
    configured_files = compilation_files(database_path, repo_root)
    if arguments.all:
      files_to_check = configured_files
    else:
      changed = changed_files(repo_root, arguments.base_ref)
      if any(path.suffix in HEADER_SUFFIXES for path in changed):
        print("Header change detected; checking all configured translation units.")
        files_to_check = configured_files
      else:
        changed_sources = {
          project_relative_path(path.resolve(), repo_root)
            for path in changed
            if path.suffix in SOURCE_SUFFIXES and is_project_source(path, repo_root)
        }
        files_to_check = [
          path
          for path in configured_files
          if project_relative_path(path, repo_root) in changed_sources
        ]
  except (RuntimeError, subprocess.CalledProcessError) as error:
    print(f"Unable to select files for clang-tidy: {error}", file=sys.stderr)
    return 2

  if not files_to_check:
    print("No configured C++ translation units to check.")
    return 0

  project_roots = {
      source_root
      for path in files_to_check
      if (source_root := project_source_root(path, repo_root)) is not None
  }
  command_prefix = [
      clang_tidy,
      f"-p={build_dir}",
      "--warnings-as-errors=*",
      f"--header-filter={arguments.header_filter or default_header_filter(project_roots)}",
  ]
  if arguments.fix:
    command_prefix.append("--fix")

  failed = False
  for source_file in files_to_check:
    command = [*command_prefix, str(source_file)]
    print("Checking", project_relative_path(source_file, repo_root))
    if subprocess.run(command, cwd=repo_root).returncode != 0:
      failed = True

  return 1 if failed else 0


if __name__ == "__main__":
  raise SystemExit(main())