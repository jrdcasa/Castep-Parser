import os
import fnmatch


# =============================================================================
def find_files_matching_pattern(root_dir_list, pattern):
    matched_files = []

    for idir in root_dir_list:
        # Walk through all directories and subdirectories recursively
        for root, dirs, files in os.walk(idir):
            # Check each file in the current directory for a match against the pattern
            for filename in files:
                if fnmatch.fnmatch(filename, pattern):
                    # If the filename matches the pattern, add its full path to the list
                    matched_files.append(os.path.join(root, filename))

    return matched_files
