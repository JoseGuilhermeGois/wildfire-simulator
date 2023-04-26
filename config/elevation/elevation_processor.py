from typing import TextIO
import math
import numpy as np

from config.business_processor import BusinessConfigProcessor, skip_lines
from .landscape import Landscape, Shape, Location
from .landscape_processor import LandscapeProcessor
from .elevation import Elevation, Shape, Location


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

        slope_steepness = self.calculate_slope_steepness()
        upslope_direction = self.calculate_upslope()

        return Elevation(slope_steepness, upslope_direction)

    @staticmethod
    def get_parameter_value(file: TextIO, definition: str) -> str:
        line = file.readline().rstrip()
        var_name, var_value = line.split()

        if var_name != definition:
            raise NameError(f"{definition} is missing!")

        return var_value

    # Unit: fraction
    def zex(self) -> float:
        """horizontal value for slope steepness and aspect calculation"""
        return -(self.elevations_horizontal_neighbors[1] - self.elevations_horizontal_neighbors[0]) / \
               2 * (Landscape.element_size / 3.281)

    # Unit: fraction
    def zey(self) -> float:
        """vertical value for slope steepness and aspect calculation"""
        return -(self.elevation_vertical_neighbors[1] - self.elevations_horizontal_neighbors[0]) / \
               2 * (Landscape.element_size / 3.281)

    # Unit: radian (fraction)
    def calculate_slope_steepness(self) -> float:
        """Calculation of slope steepness from elevation"""
        zex2 = self.zex() ** 2
        zey2 = self.zey() ** 2
        if zex2 == 0 and zey2 == 0:  # TODO: check this clause later
            return 0
        else:
            return math.sqrt(zex2 + zey2)

    # Unit: degree (azimuth)
    def calculate_aspect(self) -> float:
        """Calculation of aspect from elevation"""
        if self.zex() == 0 or self.zey() == 0:  # TODO: check this clause later
            aspect = 0
        else:
            aspect = 180 - np.arctan(self.zey() / self.zex()) + 90 * (self.zex() / abs(self.zex()))

        aspect = 180 - aspect  # to fix coordination difference between numpy array and cardinal direction
        aspect += 360 if aspect < 0 else 0  # to get positive degree
        return aspect

    # Unit: degree
    def calculate_upslope_direction(self) -> float:
        """Up-slope direction is the opposite of aspect"""
        upslope_direction = self.aspect - 180
        upslope_direction += 360 if self.upslope_direction < 0 else 0
        return upslope_direction


