from fire import Element


class TerrainTopographyBuilder:
    def __init__(self, shape: Shape, terrain_topography_facade: TerrainTopographyBuilder):
        self.shape: Shape = shape
        self.terrain_facade: TerrainTopographyBuilder = terrain_topography_facade

    def build(self):



