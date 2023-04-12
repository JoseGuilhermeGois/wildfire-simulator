from typing import TextIO

from config import BusinessConfigProcessor, skip_lines


class IgnitionsProcessor(BusinessConfigProcessor):

    def process(self, file: TextIO) -> list:
        skip_lines(file)

        return file.readline().split(',')
