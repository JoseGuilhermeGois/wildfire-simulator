from typing import TextIO

from config.business_processor import BusinessConfigProcessor, skip_lines
from .landscape import Landscape, Shape, Location


WIDTH_DEFINITION = "ncols"
LENGTH_DEFINITION = "nrows"
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
            width=int(self.get_parameter_value(file, WIDTH_DEFINITION)),
            length=int(self.get_parameter_value(file, LENGTH_DEFINITION)))

        location = Location(
            real_latitude=float(self.get_parameter_value(file, LATITUDE_COORDINATE)),
            real_longitude=float(self.get_parameter_value(file, LONGITUDE_COORDINATE)))

        element_size = int(self.get_parameter_value(file, ELEMENT_SIZE_DEFINITION))

        skip_lines(file)

        fuel_model_distribution = self.get_fuel_values(file, shape.length, shape.width)

        return Landscape(shape, location, element_size, fuel_model_distribution)

    @staticmethod
    def get_parameter_value(file: TextIO, definition: str) -> str:
        line = file.readline().rstrip()
        var_name, var_value = line.split()

        if var_name != definition:
            raise NameError(f"{definition} is missing!")

        return var_value

    def get_fuel_values(self, file: TextIO, length: int, width: int):
        fuel_values = []
        for _ in range(length + 6):
            fuel_line_values = file.readline().rsplit()[:width]
            fuel_values.append(fuel_line_values)

        return fuel_values



