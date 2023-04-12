from typing import TextIO

from config import BusinessConfigProcessor, skip_lines


class IterationsProcessor(BusinessConfigProcessor):

    def process(self, file: TextIO):

        skip_lines(file)

        return int(file.readline().split()[0])
