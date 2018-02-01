import shutil
from os import listdir, remove
from os.path import isfile, join
"""
EXAMPLE USAGE:
bed_directory_path = "/Users/saureous/data/bed_alleles_copy"
output_file_path = "output/bed_output/concatenated.bed"

file_manager = FileManager()
file_paths = file_manager.get_filenames_from_directory(directory_path=bed_directory_path)
file_manager.concatenate_files(file_paths=file_paths, output_file_path=output_file_path)
file_manager.delete_files(file_paths=file_paths)
"""

class FileManager:
    """
    Does simple file operations like concatenation, fetching a list of paths for files in a directory, deleting files
    """
    def concatenate_files(self, file_paths, output_file_path):
        with open(output_file_path, 'wb') as out_file:
            for file_path in file_paths:
                with open(file_path, 'rb') as in_file:
                    # 10MB per writing chunk to avoid reading big file into memory.
                    shutil.copyfileobj(in_file, out_file, 1024*1024*10)

    def get_filenames_from_directory(self, directory_path):
        file_paths = [join(directory_path, file) for file in listdir(directory_path) if isfile(join(directory_path, file))]
        return file_paths

    def delete_files(self, file_paths):
        for file_path in file_paths:
            remove(file_path)

