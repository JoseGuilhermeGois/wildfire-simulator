from typing import TextIO

from config import BusinessConfigProcessor
from config.environment.environment_facade import LoadWindFacade


class EnvironmentProcessor(BusinessConfigProcessor):

    def __init__(self, wind_facade: LoadWindFacade):
        self.wind_facade = wind_facade

    def process(self, file: TextIO):

        self.wind_facade.get_wind(file)
