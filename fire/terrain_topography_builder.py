from fire.element import Element
from fire.terrain_topography_facade import TerrainTopographyFacade
from config.landscape import Shape


class TerrainTopographyBuilder:
    def __init__(self, shape: Shape, terrain_topography_facade: TerrainTopographyFacade):
        self.shape: Shape = shape
        self.terrain_facade: TerrainTopographyFacade = terrain_topography_facade

    def build(self) -> list[list[Element]]:
        return [self.get_longitude_line_elements(longitude) for longitude in range(self.shape.width)]

    def get_longitude_line_elements(self, longitude):
        return [self.terrain_facade.get_element(longitude, latitude) for latitude in range(self.shape.length)]
