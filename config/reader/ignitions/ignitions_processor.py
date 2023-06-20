from typing import TextIO

from config.reader.provider import BusinessConfigProcessor, skip_lines


class IgnitionsProcessor(BusinessConfigProcessor):

    def process(self, file: TextIO) -> list[list]:

        skip_lines(file)

        return [file.readline().split() for _ in file]
