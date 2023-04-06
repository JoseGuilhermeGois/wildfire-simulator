from typing import TextIO

from config.business_processor import BusinessConfigProcessor, skip_lines
from landscape import Landscape, Shape, Location


LENGTH_DEFINITION = "ncols"
WIDTH_DEFINITION = "nrows"
LATITUDE_COORDINATE = "xllcorner"
LONGITUDE_COORDINATE = "yllcorner"
ELEMENT_SIZE_DEFINITION = "cellsize"


class LandscapeProcessor(BusinessConfigProcessor):

    def process(self, file: TextIO) -> Landscape:
        """
        We use the fuel model distribution file to discover some informations of the landscape. We discover
        the width and length of the landscape, the cell size and initial coordinates.
        :return: Landscape details
        """
        shape = Shape(
            length=self.get_parameter_value(file, LENGTH_DEFINITION),
            width=self.get_parameter_value(file, WIDTH_DEFINITION))

        location = Location(
            latitude=self.get_parameter_value(file, LATITUDE_COORDINATE),
            longitude=self.get_parameter_value(file, LONGITUDE_COORDINATE))

        element_size = self.get_parameter_value(file, ELEMENT_SIZE_DEFINITION)

        skip_lines(file)

        fuel_model_distribution = [line.rsplit() for line in file]

        return Landscape(shape, location, element_size, fuel_model_distribution)

    @staticmethod
    def get_parameter_value(file: TextIO, definition: str) -> int:
        var_name, var_value = file.readline().rstrip()

        if var_name != definition:
            raise NameError(f"{definition} is missing!")

        return int(var_value)


