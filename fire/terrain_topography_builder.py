from fire import Element
from fire.terrain_topography_facade import TerrainTopographyFacade


class TerrainTopographyBuilder:
    def __init__(self, shape: Shape, terrain_topography_facade: TerrainTopographyFacade):
        self.shape: Shape = shape
        self.terrain_facade: TerrainTopographyFacade = terrain_topography_facade

    def build(self):



