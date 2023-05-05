from typing import TextIO

from config import BusinessConfigProcessor, Landscape
from config.environment.environment_facade import LoadDefaults, LoadModels
from config.environment.environment import EnvironmentVariables, WindIndex


class EnvironmentProcessor(BusinessConfigProcessor):

    def __int__(self, landscape: Landscape):
        self.landscape: Landscape = landscape

    def process(self, file: TextIO):

        wind_model_name = file.readline().strip()

        if wind_model_name[:6] == "Defaults":
            wind_defaults = LoadDefaults(self.landscape.shape)
            return wind_defaults

        elif wind_model_name[:6] == "Canyon":
            nodes = list(map(float, file.readline().split()))
            EnvironmentVariables.NiWnd, EnvironmentVariables.NjWnd, EnvironmentVariables.NkWnd = list(
                map(int, nodes[:3]))
            EnvironmentVariables.AuxSingle = nodes[3]
            LoadModels.wind_matrix_dimensioning()
            WindIndex.X_Wnd[0], WindIndex.Y_Wnd[0], WindIndex.DxWnd = map(float,
                                                                          file.readline().split(","))
            WindIndex.DyWnd = WindIndex.DxWnd
            LoadModels.remove_fictitious_nodes()

            for i in range(EnvironmentVariables.NiWnd):
                for j in range(EnvironmentVariables.NjWnd):
                    for k in range(EnvironmentVariables.NkWnd):
                        WindIndex.U[i][j][k], WindIndex.V[i][j][k], WindIndex.W[i][j][k], teijk, edijk, tijk, pijk, \
                            WindIndex.Z_Wnd[i][j][k] = map(float, file.readline().split(","))

            for i in range(1, EnvironmentVariables.NiWnd):
                WindIndex.X_Wnd[i] = WindIndex.X_Wnd[0] + i * WindIndex.DxWnd

            for j in range(1, EnvironmentVariables.NjWnd):
                WindIndex.Y_Wnd[j] = WindIndex.Y_Wnd[0] + j * WindIndex.DyWnd

        else:
            nodes = map(int, file.readline().split())
            EnvironmentVariables.NiWnd, EnvironmentVariables.NjWnd, EnvironmentVariables.NkWnd = nodes
            LoadModels.wind_matrix_dimensioning()

            LoadModels.remove_fictitious_nodes()

            for i in range(EnvironmentVariables.NiWnd):
                for j in range(EnvironmentVariables.NjWnd):
                    for k in range(EnvironmentVariables.NkWnd):
                        WindIndex.X_Wnd[i], WindIndex.Y_Wnd[j], WindIndex.Z_Wnd[i][j][k] = \
                            map(float, file.readline().split())
                        WindIndex.U[i][j][k], WindIndex.V[i][j][k], WindIndex.W[i][j][k] = \
                            map(float, file.readline().split())

            WindIndex.DxWnd = WindIndex.X_Wnd[2] - WindIndex.X_Wnd[1]
            WindIndex.DyWnd = WindIndex.Y_Wnd[2] - WindIndex.Y_Wnd[1]

        LoadModels.final_calculations()
