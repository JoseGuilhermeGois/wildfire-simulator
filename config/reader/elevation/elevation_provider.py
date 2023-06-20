from typing import TextIO

from config.reader.provider import BusinessConfigProcessor, skip_lines
from .elevation import Elevation
from .elevation_facade import ElevationFacade
from ... import Landscape

WIDTH_DEFINITION = "ncols"
LENGTH_DEFINITION = "nrows"
LATITUDE_COORDINATE = "xllcorner"
LONGITUDE_COORDINATE = "yllcorner"
ELEMENT_SIZE_DEFINITION = "cellsize"
NO_DATA_VALUE = "NODATAVALUE"


class ElevationProvider(BusinessConfigProcessor):

    def __init__(self, landscape: Landscape):
        self.landscape: Landscape = landscape

    def process(self, file: TextIO) -> Elevation:
        """
        We use the elevation file to calculate slope steepness and up-slope direction. We discover
        the width and length of the elevation, the cell size and initial coordinates to compare with the landscape file.
        :return: Elevation details
        """
        width = int(self.get_parameter_value(file, WIDTH_DEFINITION))
        length = int(self.get_parameter_value(file, LENGTH_DEFINITION))

        # latitude = float(self.get_parameter_value(file, LATITUDE_COORDINATE))
        # longitude = float(self.get_parameter_value(file, LONGITUDE_COORDINATE))

        # element_size = int(self.get_parameter_value(file, ELEMENT_SIZE_DEFINITION))
        # no_data_value = int(self.get_parameter_value(file, NO_DATA_VALUE))

        for _ in range(4):
            skip_lines(file)

        elevation_reader = [line.rsplit() for line in file]
        elevation_distribution = elevation_reader[-length:]
        elevation_values = [row[:width] for row in elevation_distribution]
        correct_elevation_values = []
        for i in range(len(elevation_values) - 1, -1, -1):
            correct_elevation_values.append(elevation_distribution[i])

        results = ElevationFacade.calc_slope(
            number_of_columns=width,
            number_of_rows=length,
            element_size=self.landscape.element_size,
            elevation_distribution=elevation_distribution
        )

        return Elevation(results[0], results[1])

    def get_parameter_value(self, file: TextIO, definition: str) -> str:
        line = file.readline().rstrip()
        var_name, var_value = line.split()

        if var_name != definition:
            raise NameError(f"{definition} is missing!")

        return var_value

