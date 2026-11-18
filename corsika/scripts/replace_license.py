import os
import re

def is_text_file(file_path):
    """Check if a file is a text file."""
    try:
        with open(file_path, 'rb') as file:
            # Read a small portion of the file to check for binary content
            chunk = file.read(1024)
        return not any(b > 127 for b in chunk)  # Non-ASCII bytes indicate binary
    except Exception as e:
        print(f"Error checking file type for {file_path}: {e}")
        return False

# Define the new license header template (without the year)
NEW_LICENSE_TEMPLATE = """/*
 * (c) Copyright {year} CORSIKA Project, corsika8@kit.edu
 *
 * This software is distributed under the terms of the 3-clause BSD license.
 * See file LICENSE for a full version of the license.
 */"""

# Define the old license header pattern (simplified for matching)
OLD_LICENSE_PATTERN = re.compile(r"""/\*\s*(?:-\*-[^\n]*-\*-\s*)?
 \* \(c\) Copyright (\d{4}) CORSIKA Project, corsika8@kit.edu
 \*\s*
 \* This software is distributed under the terms of the GNU General Public
 \* Licence version 3 \(GPL Version 3\)\. See file LICENSE for a full version of
 \* the license\.
 \*/""", re.MULTILINE)

# Blacklist configuration
BLACKLIST_DIRS = [".git", "build", "venv"]
BLACKLIST_FILES = ["README.md", "LICENSE"]


def should_skip(path, root):
    """Check if a file or directory is blacklisted."""
    # Check if the path matches any blacklisted directory
    for blacklisted_dir in BLACKLIST_DIRS:
        if os.path.commonpath([os.path.abspath(os.path.join(root, blacklisted_dir)), os.path.abspath(path)]) == os.path.abspath(os.path.join(root, blacklisted_dir)):
            return True

    # Check if the file is blacklisted
    if os.path.basename(path) in BLACKLIST_FILES:
        return True

    return False


def replace_license_in_file(file_path):
    """Replace the license header in a given file."""
    try:
        if not is_text_file(file_path):
            print(f"Skipping binary file: {file_path}")
            return

        with open(file_path, "r", encoding="utf-8") as file:
            content = file.read()

        match = OLD_LICENSE_PATTERN.search(content)
        if match:
            # Extract the year from the old license
            year = match.group(1)
            # Generate the new license header with the same year
            new_license = NEW_LICENSE_TEMPLATE.format(year=year)
            # Replace the old license with the new one
            content = content.replace(match.group(0), new_license)

            with open(file_path, "w", encoding="utf-8") as file:
                file.write(content)

            print(f"Updated license in: {file_path}")
        else:
            print(f"No matching license found in: {file_path}")

    except Exception as e:
        print(f"Error processing file {file_path}: {e}")


def process_directory(directory):
    """Recursively process a directory, replacing license headers."""
    for root, dirs, files in os.walk(directory):
        # Filter out blacklisted directories
        dirs[:] = [d for d in dirs if not should_skip(os.path.join(root, d), directory)]

        for file in files:
            file_path = os.path.join(root, file)

            # Skip blacklisted files
            if should_skip(file_path, directory):
                continue

            # Process the file to replace the license
            replace_license_in_file(file_path)


if __name__ == "__main__":
    # Set the starting directory for processing
    project_directory = os.getcwd()  # Use the current working directory

    print(f"Starting license update in directory: {project_directory}")
    process_directory(project_directory)
    print("License update completed.")

