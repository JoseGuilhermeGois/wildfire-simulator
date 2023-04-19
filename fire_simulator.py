"""This file is to run the simulation."""

from config import BusinessConfigProcessor, LandscapeProcessor, FuelModelsProcessor, BaseFuelFacade, DefaultsFuel, \
    IterationsProcessor, IgnitionsProcessor, EnvironmentProcessor, LoadDefaults
from fire import State, BaseTerrainTopographyFacade, TerrainTopographyBuilder, Fire


class FireSimulator:

    def __init__(self,
                 landscape_processor: BusinessConfigProcessor,
                 fuel_models_processor: BusinessConfigProcessor,
                 environment_processor: BusinessConfigProcessor,
                 ignition_processor: BusinessConfigProcessor,
                 iterations_processor: BusinessConfigProcessor):

        self.landscape_processor = landscape_processor
        self.fuel_models_processor = fuel_models_processor
        self.environment_processor = environment_processor
        self.ignition_processor = ignition_processor
        self.iterations_processor = iterations_processor

    def run(self, landscape_filename: str, fuel_models_filename: str, environment_filename: str,
            ignition_points_filename: str, iterations_filename: str):

        landscape = self.landscape_processor.read(landscape_filename)
        fuel_models = self.fuel_models_processor.read(fuel_models_filename)
        environment = self.environment_processor.read(environment_filename)
        ignition_points = self.ignition_processor.read(ignition_points_filename)
        iterations = self.iterations_processor.read(iterations_filename)

        terrain_topography_facade = BaseTerrainTopographyFacade(landscape, fuel_models, environment)
        terrain_topography = TerrainTopographyBuilder(landscape.shape, terrain_topography_facade).build()

        for longitude, latitude in ignition_points:
            terrain_topography[longitude][latitude].state = State.BURNING

        Fire(landscape, terrain_topography, iterations).start()


if __name__ == '__main__':
    defaults_fuel = DefaultsFuel()

    fire_simulator = FireSimulator(
        landscape_processor=LandscapeProcessor(),
        fuel_models_processor=FuelModelsProcessor(BaseFuelFacade(defaults_fuel)),
        environment_processor=EnvironmentProcessor(LoadDefaults(LandscapeProcessor())),
        ignition_processor=IgnitionsProcessor(),
        iterations_processor=IterationsProcessor()
    )

    fire_simulator.run(
        landscape_filename="resources/fuel_1m_gestosa.asc",
        fuel_models_filename="resources/fuelmodelsTest.fls",
        environment_filename="resources/default_wind.asc",
        ignition_points_filename="resources/ignitions.asc",
        iterations_filename="resources/time_step.asc"
    )
