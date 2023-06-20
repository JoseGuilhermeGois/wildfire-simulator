from typing import TextIO

from config.reader.provider import BusinessConfigProcessor
from config.reader.landscape import Landscape
from config.reader.environment.environment_facade import LoadModels
from config.reader.environment.environment import WindIndex, EnvironmentValues


class EnvironmentProvider(BusinessConfigProcessor):

    def __init__(self, landscape: Landscape):
        self.landscape: Landscape = landscape

    def process(self, file: TextIO) -> WindIndex:

        elevation = [[EnvironmentValues(slope_steepness=20, aspect=300) for _ in range(self.landscape.shape.width)] for _ in range(
            self.landscape.shape.length)]

        wind_model_name = file.readline().strip()
        wind_index = WindIndex

        if wind_model_name[:6] == "Canyon":
            nodes = list(map(float, file.readline().split()))
            wind_index.NiWnd, wind_index.NjWnd, wind_index.NkWnd = list(map(int, nodes[:3]))
            wind_index.AuxSingle = nodes[3]
            LoadModels.wind_matrix_dimensioning(wind_index)
            wind_index.X_Wnd[0], wind_index.Y_Wnd[0], wind_index.DxWnd = list(map(float,
                                                                              file.readline().split(",")))
            wind_index.DyWnd = wind_index.DxWnd
            LoadModels.remove_fictitious_nodes(wind_index)

            for i in range(int(wind_index.NiWnd)):
                for j in range(int(wind_index.NjWnd)):
                    for k in range(int(wind_index.NkWnd)):
                        wind_index.U[i][j][k], wind_index.V[i][j][k], wind_index.W[i][j][k], teijk, edijk, tijk, pijk, \
                            wind_index.Z_Wnd[i][j][k] = list(map(float, file.readline().split(",")))

            for i in range(1, wind_index.NiWnd):
                wind_index.X_Wnd[i] = wind_index.X_Wnd[0] + i * wind_index.DxWnd

            for j in range(1, wind_index.NjWnd):
                wind_index.Y_Wnd[j] = wind_index.Y_Wnd[0] + j * wind_index.DyWnd

        else:
            nodes = map(int, file.readline().split())
            wind_index.NiWnd, wind_index.NjWnd, wind_index.NkWnd = nodes
            LoadModels.wind_matrix_dimensioning(wind_index)

            LoadModels.remove_fictitious_nodes(wind_index)

            for i in range(wind_index.NiWnd):
                for j in range(wind_index.NjWnd):
                    for k in range(wind_index.NkWnd):
                        wind_index.X_Wnd[i], wind_index.Y_Wnd[j], wind_index.Z_Wnd[i][j][k] = \
                            map(float, file.readline().split())
                        wind_index.U[i][j][k], wind_index.V[i][j][k], wind_index.W[i][j][k] = \
                            map(float, file.readline().split())

            wind_index.DxWnd = wind_index.X_Wnd[2] - wind_index.X_Wnd[1]
            wind_index.DyWnd = wind_index.Y_Wnd[2] - wind_index.Y_Wnd[1]

        LoadModels.final_calculations(wind_index)

        return WindIndex(
            wind_model_name,
            wind_index.U,
            wind_index.V,
            wind_index.W,
            wind_index.Z_Wnd,
            wind_index.ZZ_Wnd,
            wind_index.X_Wnd,
            wind_index.XX_Wnd,
            wind_index.Y_Wnd,
            wind_index.YY_Wnd,
            wind_index.DxWnd,
            wind_index.DyWnd,
            wind_index.X_IniWnd,
            wind_index.Y_IniWnd,
            wind_index.X_FinWnd,
            wind_index.Y_FinWnd,
            wind_index.NiWnd,
            wind_index.NjWnd,
            wind_index.NkWnd,
            wind_index.AuxSingle)
