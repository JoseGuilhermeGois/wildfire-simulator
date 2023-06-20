from typing import Protocol
from typing import TextIO, TypeVar


T = TypeVar("T")


class BusinessConfigProcessor(Protocol):

    def read(self, filename: str) -> T:
        try:
            file = open(file=filename, mode='r', encoding='latin-1')
        except FileNotFoundError:
            raise f"File {filename} not found. Aborting!"
        except OSError:
            raise f"OS error occurred trying to open {filename}!"
        except Exception as err:
            raise (f"Unexpected error opening {filename} is", repr(err))
        else:
            with file:
                return self.process(file)

    def process(self, file: TextIO) -> T:
        ...


def skip_lines(file: TextIO, lines: int = 1):
    for _ in range(lines):
        next(file)

