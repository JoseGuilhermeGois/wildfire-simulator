"""This file is to run the simulation."""


class FireSimulator:

    def __init__(self, landscape_processor, fuel_models_processor):
        self.landscape_processor = landscape_processor
        self.fuel_models_processor = fuel_models_processor

    def start(self, landscape_filename: str, fuel_models_filename: str, ignition_points: tuple) -> str:
        landscape = self.landscape_processor.read(landscape_filename)
        fuel_models_filename = self.fuel_models_processor.read(fuel_models_filename)

