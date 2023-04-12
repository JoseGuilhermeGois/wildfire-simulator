from typing import TextIO

from config import BusinessConfigProcessor
from config.environment.environment_facade import LoadCanyon, LoadNuatmos
from config.landscape import Landscape


class LoadWind(BusinessConfigProcessor):

    def process(self, file: TextIO):

        wind_model_name = file.readline().strip()

        if wind_model_name[:6] == "Canyon":
            LoadCanyon()
        else:
            LoadNuatmos()


def create_fake_environment(landscape: Landscape, default_environment: Environment):

    return[[default_environment for i in range(landscape.shape.width) for i in range(landscape.shape.length)]]
