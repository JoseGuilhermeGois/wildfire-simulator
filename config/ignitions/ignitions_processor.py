from typing import TextIO

from config import BusinessConfigProcessor, skip_lines


class IgnitionsProcessor(BusinessConfigProcessor):

    def process(self, file: TextIO) -> list:

        return [file.readline().split(',') for 1 in file]
