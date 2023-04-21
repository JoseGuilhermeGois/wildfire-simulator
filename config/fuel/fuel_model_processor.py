from typing import TextIO

from config.business_processor import BusinessConfigProcessor, skip_lines
from .fuel import Fuel
from .fuel_facade import FuelFacade


BATCH_OF_LINES_FOR_FUEL_MODEL_CHARACTERISTICS = 19


class FuelModelsProcessor(BusinessConfigProcessor):

    def __init__(self, fuel_facade: FuelFacade):
        self.fuel_facade: FuelFacade = fuel_facade

    def process(self, file: TextIO) -> dict[str, Fuel]:

        skip_lines(file)

        fuel_models = {}
        counter = 0
        while counter < 8:
            fuel_model_name = self.read_fuel_model_metadata(file, counter)
            fuel_model_characteristics = self.fuel_model_characteristics(file)
            fuel_model_id = str(counter)

            if len(fuel_model_characteristics) < BATCH_OF_LINES_FOR_FUEL_MODEL_CHARACTERISTICS:
                break

            fuel = self.fuel_facade.get_fuel_models(fuel_model_name, fuel_model_characteristics)
            fuel_models[fuel_model_id] = fuel

            counter += 1

        return fuel_models

    @staticmethod
    def read_fuel_model_metadata(file: TextIO, counter: int) -> str | None:

        if counter == 0:
            fuel_model_metadata = file.readline().split("-")[1: 3]
        else:
            fuel_model_metadata = file.readline().split()

        if not fuel_model_metadata:
            return None

        fuel_model_name = " ".join(fuel_model_metadata)
        return fuel_model_name

    @staticmethod
    def fuel_model_characteristics(file: TextIO) -> list[float]:
        return [float(file.readline().split()[0]) for _ in range(BATCH_OF_LINES_FOR_FUEL_MODEL_CHARACTERISTICS)]
