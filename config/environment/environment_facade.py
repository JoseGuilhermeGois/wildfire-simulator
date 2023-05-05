import math
from dataclasses import dataclass
from typing import TextIO, Protocol

import numpy as np

from config.environment.environment import WindIndex, Environment, EnvironmentVariables
from config.business_processor import skip_lines
from config.landscape import Shape


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

    @staticmethod
    def wind_matrix_dimensioning():
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

    @staticmethod
    def remove_fictitious_nodes():
        # The index 0 corresponds to the fictitious nodes.
        # Thus, the mesh defined initially corresponds, in i, from 1 to NiWnd
        EnvironmentVariables.NiWnd -= 1
        EnvironmentVariables.NjWnd -= 1
        EnvironmentVariables.NkWnd -= 1

    @staticmethod
    def final_calculations():
        for i in range(EnvironmentVariables.NiWnd - 1):
            WindIndex.XX_Wnd[i] = 0.5 * (WindIndex.X_Wnd[i] + WindIndex.X_Wnd[i + 1])

        for j in range(EnvironmentVariables.NjWnd - 1):
            WindIndex.YY_Wnd[j] = 0.5 * (WindIndex.Y_Wnd[j] + WindIndex.Y_Wnd[j + 1])

        WindIndex.X_IniWnd = WindIndex.X_Wnd[1]
        WindIndex.Y_IniWnd = WindIndex.Y_Wnd[1]
        WindIndex.X_FinWnd = WindIndex.X_Wnd[EnvironmentVariables.NiWnd]
        WindIndex.Y_FinWnd = WindIndex.Y_Wnd[EnvironmentVariables.NjWnd]

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
