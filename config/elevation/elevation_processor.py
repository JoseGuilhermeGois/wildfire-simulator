from typing import TextIO
import math
import numpy as np

from config.business_processor import BusinessConfigProcessor, skip_lines
from .landscape import Landscape, Shape, Location
from .landscape_processor import LandscapeProcessor
from .elevation import Elevation


WIDTH_DEFINITION = "ncols"
LENGTH_DEFINITION = "nrows"
LATITUDE_COORDINATE = "xllcorner"
LONGITUDE_COORDINATE = "yllcorner"
ELEMENT_SIZE_DEFINITION = "cellsize"


class ElevationProcessor(BusinessConfigProcessor):

    def __init__(self, landscape: Landscape):
        self.landscape = landscape

    def process(self, file: TextIO) -> Elevation:
        """
        We use the elevation file to calculate slope steepness and up-slope direction. We discover
        the width and length of the elevation, the cell size and initial coordinates to compare with the landscape file.
        :return: Elevation details
        """
        width = int(self.get_parameter_value(file, WIDTH_DEFINITION))
        length = int(self.get_parameter_value(file, LENGTH_DEFINITION))

        latitude = float(self.get_parameter_value(file, LATITUDE_COORDINATE))
        longitude = float(self.get_parameter_value(file, LONGITUDE_COORDINATE))

        element_size = int(self.get_parameter_value(file, ELEMENT_SIZE_DEFINITION))

        # if width == landscape.shape.width:
        #     continue
        # elif length == landsape.shape.length:
        #
        #      raise ValueError


        skip_lines(file)

        elevation_reader = [line.rsplit() for line in file]
        elevation_distribution = elevation_reader[:length]

        for i in range(len(elevation_distribution)):
            elevation_distribution[i] = elevation_distribution[i][:width]

        return Elevation(elevation_distribution)

    @staticmethod
    def get_parameter_value(file: TextIO, definition: str) -> str:
        line = file.readline().rstrip()
        var_name, var_value = line.split()

        if var_name != definition:
            raise NameError(f"{definition} is missing!")

        return var_value


