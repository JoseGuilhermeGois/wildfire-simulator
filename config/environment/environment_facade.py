import math
from dataclasses import dataclass
from typing import TextIO, Protocol

import numpy as np

from config.environment.environment import WindIndex, Environment
from config.business_processor import skip_lines
from config.landscape import Shape


@dataclass
class EnvironmentVariables:
    NiWnd: int
    NjWnd: int
    NkWnd: int
    AuxSingle: float
    DxWnd: float
    DyWnd: float
    X_IniWnd: float
    Y_IniWnd: float
    X_FinWnd: float
    Y_FinWnd: float


class LoadWindFacade(Protocol):

    def get_wind(self, file: TextIO) -> list:
        ...


class LoadDefaults(LoadWindFacade):

    def __init__(self, shape: Shape):
        self.shape: Shape = shape

    def get_wind(self, file: TextIO) -> list[list[Environment]]:

        skip_lines(file)

        wind_speed = round(float(file.readline().split()[0]), 4) * 196.9
        wind_direction = 180 - round(float(file.readline().split()[0]), 4)
        slope_steepness = math.tan(math.radians(round(float(file.readline().split()[0]), 4)))
        aspect = 180 - round(float(file.readline().split()[0]), 4)

        default_environment = Environment(wind_speed=wind_speed,
                                          wind_direction=wind_direction,
                                          slope_steepness=slope_steepness,
                                          aspect=aspect)

        return [[default_environment for _ in range(self.shape.width)] for _ in range(
            self.shape.length)]


class LoadModels(LoadWindFacade):

    def get_wind(self, file: TextIO):

        wind_model_name = file.readline().strip()

        if wind_model_name[:6] == "Canyon":
            nodes = list(map(float, file.readline().split()))
            EnvironmentVariables.NiWnd, EnvironmentVariables.NjWnd, EnvironmentVariables.NkWnd = list(map(int, nodes[:3]))
            EnvironmentVariables.AuxSingle = nodes[3]
            self.wind_matrix_dimensioning()
            WindIndex.X_Wnd[0], WindIndex.Y_Wnd[0], EnvironmentVariables.DxWnd = map(float, file.readline().split(","))
            EnvironmentVariables.DyWnd = EnvironmentVariables.DxWnd
            self.remove_fictitious_nodes()

            for i in range(EnvironmentVariables.NiWnd):
                for j in range(EnvironmentVariables.NjWnd):
                    for k in range(EnvironmentVariables.NkWnd):
                        WindIndex.U[i][j][k], WindIndex.V[i][j][k], WindIndex.W[i][j][k], teijk, edijk, tijk, pijk, \
                            WindIndex.Z_Wnd[i][j][k] = map(float, file.readline().split(","))

            for i in range(1, EnvironmentVariables.NiWnd):
                WindIndex.X_Wnd[i] = WindIndex.X_Wnd[0] + i * EnvironmentVariables.DxWnd

            for j in range(1, EnvironmentVariables.NjWnd):
                WindIndex.Y_Wnd[j] = WindIndex.Y_Wnd[0] + j * EnvironmentVariables.DyWnd

        else:
            nodes = map(int, file.readline().split())
            EnvironmentVariables.NiWnd, EnvironmentVariables.NjWnd, EnvironmentVariables.NkWnd = nodes
            self.wind_matrix_dimensioning()

            self.remove_fictitious_nodes()

            for i in range(EnvironmentVariables.NiWnd):
                for j in range(EnvironmentVariables.NjWnd):
                    for k in range(EnvironmentVariables.NkWnd):
                        WindIndex.X_Wnd[i], WindIndex.Y_Wnd[j], WindIndex.Z_Wnd[i][j][k] = \
                            map(float, file.readline().split())
                        WindIndex.U[i][j][k], WindIndex.V[i][j][k], WindIndex.W[i][j][k] = \
                            map(float, file.readline().split())

            EnvironmentVariables.DxWnd = WindIndex.X_Wnd[2] - WindIndex.X_Wnd[1]
            EnvironmentVariables.DyWnd = WindIndex.Y_Wnd[2] - WindIndex.Y_Wnd[1]

        results = self.final_calculations()

        return results

    def wind_matrix_dimensioning(self):
        # Create and initialize arrays
        WindIndex.X_Wnd = np.zeros(EnvironmentVariables.NiWnd + 1)
        WindIndex.Y_Wnd = np.zeros(EnvironmentVariables.NjWnd + 1)
        WindIndex.Z_Wnd = np.zeros((EnvironmentVariables.NiWnd + 1, EnvironmentVariables.NjWnd + 1,
                                    EnvironmentVariables.NkWnd + 1))

        WindIndex.XX_Wnd = np.zeros(EnvironmentVariables.NiWnd)
        WindIndex.YY_Wnd = np.zeros(EnvironmentVariables.NjWnd)
        WindIndex.ZZ_Wnd = np.zeros((EnvironmentVariables.NiWnd + 1, EnvironmentVariables.NjWnd + 1,
                                     EnvironmentVariables.NkWnd + 1))

        WindIndex.U = np.zeros((EnvironmentVariables.NiWnd + 1, EnvironmentVariables.NjWnd + 1,
                                EnvironmentVariables.NkWnd + 1))
        WindIndex.Y = np.zeros((EnvironmentVariables.NiWnd + 1, EnvironmentVariables.NjWnd + 1,
                                EnvironmentVariables.NkWnd + 1))
        WindIndex.W = np.zeros((EnvironmentVariables.NiWnd + 1, EnvironmentVariables.NjWnd + 1,
                                EnvironmentVariables.NkWnd + 1))

    def remove_fictitious_nodes(self):
        # The index 0 corresponds to the fictitious nodes.
        # Thus, the mesh defined initially corresponds, in i, from 1 to NiWnd
        EnvironmentVariables.NiWnd -= 1
        EnvironmentVariables.NjWnd -= 1
        EnvironmentVariables.NkWnd -= 1

    def final_calculations(self):
        for i in range(EnvironmentVariables.NiWnd - 1):
            WindIndex.XX_Wnd[i] = 0.5 * (WindIndex.X_Wnd[i] + WindIndex.X_Wnd[i + 1])

        for j in range(EnvironmentVariables.NjWnd - 1):
            WindIndex.YY_Wnd[j] = 0.5 * (WindIndex.Y_Wnd[j] + WindIndex.Y_Wnd[j + 1])

        Environment.X_IniWnd = WindIndex.X_Wnd[1]
        Environment.Y_IniWnd = WindIndex.Y_Wnd[1]
        Environment.X_FinWnd = WindIndex.X_Wnd[EnvironmentVariables.NiWnd]
        Environment.Y_FinWnd = WindIndex.Y_Wnd[EnvironmentVariables.NjWnd]

        for i in range(EnvironmentVariables.NiWnd - 1):
            for j in range(EnvironmentVariables.NjWnd - 1):
                for k in range(EnvironmentVariables.NkWnd - 1):
                    WindIndex.ZZ_Wnd[i][j][k] = 0.125 * (WindIndex.Z_Wnd[i][j][k] + WindIndex.Z_Wnd[i + 1][j][k] +
                                                         WindIndex.Z_Wnd[i + 1][j + 1][k] + WindIndex.Z_Wnd[i][j + 1][
                                                             k] +
                                                         WindIndex.Z_Wnd[i][j][k + 1] + WindIndex.Z_Wnd[i + 1][j][
                                                             k + 1] +
                                                         WindIndex.Z_Wnd[i + 1][j + 1][k + 1] +
                                                         WindIndex.Z_Wnd[i][j + 1][k + 1] - 2 *
                                                         WindIndex.Z_Wnd[i][j][1] - 2 * WindIndex.Z_Wnd[i + 1][j][
                                                             1] - 2 *
                                                         WindIndex.Z_Wnd[i + 1][j + 1][k] - 2 *
                                                         WindIndex.Z_Wnd[i][j + 1][1])

        return WindIndex.XX_Wnd, WindIndex.YY_Wnd, WindIndex.ZZ_Wnd
