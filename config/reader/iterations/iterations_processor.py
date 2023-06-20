from typing import TextIO

from config.reader import BusinessConfigProcessor, skip_lines


class IterationsProcessor(BusinessConfigProcessor):

    def process(self, file: TextIO):

        skip_lines(file)

        return file.readline().split()
