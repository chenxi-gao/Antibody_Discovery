import os
import sys


def ensure_directory_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)


class Logger:
    def __init__(self, base_file_name, timestamp, mode, program_name):
        # Ensure all parameters are provided
        if not base_file_name:
            raise ValueError("Base file name cannot be null")
        if not timestamp:
            raise ValueError("Timestamp cannot be null")
        if not mode:
            raise ValueError("Mode cannot be null")
        if not program_name:
            raise ValueError("Program name cannot be null")

        current_dir = os.getcwd()
        logs_dir_path = os.path.join(current_dir, "test/logs" if mode == "test" else "logs", program_name)
        log_file_name = f"{'[test]' if mode == 'test' else ''}[{timestamp}]{base_file_name}.log"
        err_file_name = f"{'[test]' if mode == 'test' else ''}[{timestamp}]{base_file_name}.err"

        ensure_directory_exists(logs_dir_path)

        log_file_path = os.path.join(logs_dir_path, log_file_name)
        err_file_path = os.path.join(logs_dir_path, err_file_name)

        self.log_writer = open(log_file_path, 'a')
        self.err_writer = open(err_file_path, 'a')

    def log(self, message):
        print(message, file=self.log_writer)
        print(message)

    def blank(self):
        print(file=self.log_writer)
        print()

    def error(self, message):
        print(f"ERROR: {message}", file=self.err_writer)
        print(f"ERROR: {message}", file=sys.stderr)

    def close(self):
        if self.log_writer:
            self.log_writer.close()
        if self.err_writer:
            self.err_writer.close()
