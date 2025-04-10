from abc import ABC
import csv
from io import TextIOWrapper
from typing import Any, Callable

class FileProcessError(Exception, ABC):
    """Abatract class for file process error"""
    
    def __init__(self, file_path: str, message: str):
        self.file_path = file_path
        super().__init__(f"File Error [{self.file_path}]: {message}")

class FileOpenError(FileProcessError):
    """File open error"""
    
    def __init__(self, file_path: str):
        super().__init__(file_path, "Can not open file.")

class FileIOError(FileProcessError):
    """File IO error"""
    
    def __init__(self, file_path: str):
        super().__init__(file_path, "File IO error.")

def _process_file(file_path: str, mode: str, call_back: Callable[[TextIOWrapper], None]) -> Any:
    try:
        with open(file_path, mode) as f:
            return call_back(f)
    except FileNotFoundError:
        raise FileOpenError(file_path)
    except IOError:
        raise FileIOError(file_path)

def _read_file(file_path: str, call_back: Callable[[TextIOWrapper], None]) -> Any:
    return _process_file(file_path, 'r', call_back)

def read_file(file_path: str) -> tuple[str, ...]:
    def _read(f: TextIOWrapper):
        return tuple(
            line.strip() for line in f.readlines()
        )
    return _read_file(file_path, _read)

def _read_delimited_file(file_path: str, delimiter: str) -> [list[list[str]]]:
    def _read(f: TextIOWrapper):
        reader = csv.reader(f, delimiter=delimiter, lineterminator='\n')
        lines = [line for line in reader]
        return lines
    return _read_file(file_path, _read)

def read_tsv_file(file_path: str) -> [list[list[str]]]:
    return _read_delimited_file(file_path, '\t')

def _write_file(file_path: str, call_back: Callable[[TextIOWrapper], None]) -> None:
    _process_file(file_path, 'w', call_back)

def write_file(file_path: str, data: list[str]) -> None:
    def _write(f: TextIOWrapper):
        f.write("\n".join(data))
    _write_file(file_path, _write)
    print(f"File written: {file_path}")

def _write_delimited_file(file_path: str, data: list[list[str]], delimiter: str) -> None:
    def _write(f: TextIOWrapper):
        writer = csv.writer(f, delimiter=delimiter, lineterminator='\n')
        writer.writerows(data)
    _write_file(file_path, _write)

def write_tsv_file(file_path: str, data: list[list[str]]) -> None:
    _write_delimited_file(file_path, data, '\t')
