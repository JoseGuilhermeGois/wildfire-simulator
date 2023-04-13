from typing import TextIO

from config import BusinessConfigProcessor
from config.environment.environment_facade import LoadWindFacade
from config import Landscape
from environment import Environment


class EnvironmentProcessor(BusinessConfigProcessor):

    def __init__(self, wind_facade: LoadWindFacade):
        self.wind_facade = wind_facade

    def process(self, file: TextIO):

        self.wind_facade.get_wind(file)


def create_fake_environment(landscape: Landscape, default_environment: Environment):

    return[[default_environment for i in range(landscape.shape.width) for i in range(landscape.shape.length)]]
