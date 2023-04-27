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
        iterations = int(self.iterations_processor.read(iterations_filename)[0])

        terrain_topography_facade = BaseTerrainTopographyFacade(landscape, fuel_models, environment)
        terrain_topography = TerrainTopographyBuilder(landscape.shape, terrain_topography_facade).build()

        ignitions_counter = 0
        for pair in ignition_points:
            try:
                latitude, longitude = pair
            except ValueError:
                ignition_values = len(pair)
                raise RuntimeError("Length of the pair ({}) is not 2.".format(ignition_values))

            ignitions_counter += 1
            terrain_topography[int(longitude)][int(latitude)].state = State.BURNING

        Fire(landscape, terrain_topography, iterations, ignitions_counter).start()


if __name__ == '__main__':
    defaults_fuel = DefaultsFuel()

    landscape_variables = LandscapeProcessor()
    landscape_env = landscape_variables.read("resources/fuel_1m_gestosa.asc")

    fire_simulator = FireSimulator(
        landscape_processor=LandscapeProcessor(),
        fuel_models_processor=FuelModelsProcessor(BaseFuelFacade(defaults_fuel)),
        environment_processor=EnvironmentProcessor(LoadDefaults(landscape_env.shape)),
        ignition_processor=IgnitionsProcessor(),
        iterations_processor=IterationsProcessor()
    )

    fire_simulator.run(
        landscape_filename="resources/fuel_1m_gestosa.asc",
        fuel_models_filename="resources/LopesFuelModels.fls",
        environment_filename="resources/default_wind.asc",
        ignition_points_filename="resources/ignitions.asc",
        iterations_filename="resources/time_step.asc"
    )
