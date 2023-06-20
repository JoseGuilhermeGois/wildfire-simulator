"""This file is to run the simulation."""
from enum import Enum

from config.reader import BusinessConfigProcessor, LandscapeProvider, FuelModelsProcessor, BaseFuelFacade, DefaultsFuel, \
    IterationsProcessor, IgnitionsProcessor, EnvironmentProvider
from config.reader.elevation.elevation_provider import ElevationProvider
from fire import State, BaseTerrainTopographyFacade, TerrainTopographyBuilder, Fire, CombustibleElement


class ConfigFilenames(Enum):
    LANDSCAPE_FILENAME = "resources/LopesFuelsMap.asc"
    FUEL_MODELS_FILENAME = "resources/LopesFuelModels.fls"
    ENVIRONMENT_FILENAME = "resources/LopesCanyon.out"
    IGNITION_POINTS_FILENAME = "resources/ignitions.asc"
    ITERATIONS_FILENAME = "resources/time_step.asc"
    ELEVATION_FILENAME = "resources/LopesDEM.asc"


class FireSimulator:

    def __init__(self,
                 landscape_provider: BusinessConfigProcessor,
                 fuel_models_provider: BusinessConfigProcessor,
                 environment_provider: BusinessConfigProcessor,
                 ignition_provider: BusinessConfigProcessor,
                 iterations_provider: BusinessConfigProcessor,
                 elevation_provider: BusinessConfigProcessor):

        self.landscape_provider = landscape_provider
        self.fuel_models_provider = fuel_models_provider
        self.environment_provider = environment_provider
        self.ignition_provider = ignition_provider
        self.iterations_provider = iterations_provider
        self.elevation_provider = elevation_provider

        self.ignitions_list: list[CombustibleElement] = []

    def run(self, landscape_filename: str, fuel_models_filename: str, environment_filename: str,
            ignition_points_filename: str, iterations_filename: str, elevation_filename: str):

        landscape = self.landscape_provider.read(landscape_filename)
        fuel_models = self.fuel_models_provider.read(fuel_models_filename)
        environment = self.environment_provider.read(environment_filename)
        elevation = self.elevation_provider.read(elevation_filename)
        ignition_points = self.ignition_provider.read(ignition_points_filename)
        iterations = int(self.iterations_provider.read(iterations_filename)[0])

        terrain_topography_facade = BaseTerrainTopographyFacade(landscape, fuel_models, environment, elevation)
        terrain_topography = TerrainTopographyBuilder(landscape.shape, terrain_topography_facade).build()

        ignitions_counter = 0
        for pair in ignition_points:
            try:
                latitude, longitude = pair
            except ValueError:
                ignition_values = len(pair)
                raise RuntimeError("Length of the pair ({}) is not 2.".format(ignition_values))

            element = terrain_topography[int(longitude)][int(latitude)]

            if isinstance(element, CombustibleElement):
                element.state = State.BURNING
                ignitions_counter += 1
                element.time_of_ignition = 0
                self.ignitions_list.append(element)

        Fire(landscape, terrain_topography, iterations, ignitions_counter, self.ignitions_list).start()


if __name__ == '__main__':
    defaults_fuel = DefaultsFuel()

    landscape_instance = LandscapeProvider()
    environment_instance = landscape_instance.read(ConfigFilenames.LANDSCAPE_FILENAME.value)

    fire_simulator = FireSimulator(
        landscape_provider=LandscapeProvider(),
        fuel_models_provider=FuelModelsProcessor(BaseFuelFacade(defaults_fuel)),
        environment_provider=EnvironmentProvider(environment_instance),
        ignition_provider=IgnitionsProcessor(),
        iterations_provider=IterationsProcessor(),
        elevation_provider=ElevationProvider(environment_instance)
    )

    fire_simulator.run(
        landscape_filename=ConfigFilenames.LANDSCAPE_FILENAME.value,
        fuel_models_filename=ConfigFilenames.FUEL_MODELS_FILENAME.value,
        environment_filename=ConfigFilenames.ENVIRONMENT_FILENAME.value,
        ignition_points_filename=ConfigFilenames.IGNITION_POINTS_FILENAME.value,
        iterations_filename=ConfigFilenames.ITERATIONS_FILENAME.value,
        elevation_filename=ConfigFilenames.ELEVATION_FILENAME.value
    )
