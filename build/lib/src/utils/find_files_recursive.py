import os
import fnmatch


# =============================================================================
def find_files_matching_pattern(root_dir, pattern):
    matched_files = []

    # Walk through all directories and subdirectories recursively
    for root, dirs, files in os.walk(root_dir):
        # Check each file in the current directory for a match against the pattern
        for filename in files:
            if fnmatch.fnmatch(filename, pattern):
                # If the filename matches the pattern, add its full path to the list
                matched_files.append(os.path.join(root, filename))

    return matched_files
