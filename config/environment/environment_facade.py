
from dataclasses import dataclass
from typing import TextIO, Protocol

import numpy as np

from config.environment.environment import WindIndex, Environment
from config.business_processor import skip_lines
from config.landscape import Landscape


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

    def __init__(self, terrain_conditions):
        self.terrain_conditions = terrain_conditions

    def get_wind(self, file: TextIO):

        skip_lines(file)

        wind_speed = round(float(file.readline().strip()[0]), 3)
        wind_direction = round(float(file.readline().strip()[0]), 3)
        slope_steepness = round(float(file.readline().strip()[0]), 3)
        aspect = round(float(file.readline().strip()[0]), 3)

        default_environment = Environment(wind_speed=wind_speed,
                                          wind_direction=wind_direction,
                                          slope_steepness=slope_steepness,
                                          aspect=aspect)

        return self.create_fake_environment(self.terrain_conditions, default_environment)

    def create_fake_environment(self, landscape: Landscape, default_environment: Environment):
        return [[default_environment for i in range(landscape.shape.length) for i in range(landscape.shape.width)]]


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

    def interpol_wind_speed(self):
        CrdX = X_IniDtm + (i - 1) * DxDtm + 0.5 * DxDtm
        CrdY = Y_IniDtm + (j - 1) * DyDtm + 0.5 * DyDtm

        iv = int((CrdX - X_IniWnd) / DxWnd + 0.5)
        jv = int((CrdY - Y_IniWnd) / DyWnd + 0.5)

        if (iv < 0 or iv >= NiWnd - 0 or jv < 0 or jv >= NjWnd - 0):
            Vel_U, Vel_V, Vel_W = 0, 0, 0
            return

        xv1 = XX_Wnd[iv]
        xv2 = XX_Wnd[iv + 1]
        yv1 = YY_Wnd[jv]
        yv2 = YY_Wnd[jv + 1]

        Altura = ZZ_Wnd[iv, jv, 1]
        flngt = FlmLngt(ModeloComb[i, j])
        veglngt = deltaD(ModeloComb[i, j])

        # ------ Calculation of the U wind speed -------
        velfga = interpol(yv1, U[iv, jv, 1], yv2, U[iv, jv + 1, 1], CrdY)
        velfgb = interpol(yv1, U[iv + 1, jv, 1], yv2, U[iv + 1, jv + 1, 1], CrdY)
        result = interpol(xv1, velfga, xv2, velfgb, CrdX)
        Vel_U = result * (1 + 0.36 * veglngt / flngt) / (log((Altura + 0.36 * veglngt) / (0.13 * veglngt))) * (
                log((flngt / veglngt + 0.36) / 0.13) - 1)

        # ------ Calculation of the V wind speed -------
        velfga = interpol(yv1, V[iv, jv, 1], yv2, V[iv, jv + 1, 1], CrdY)
        velfgb = interpol(yv1, V[iv + 1, jv, 1], yv2, V[iv + 1, jv + 1, 1], CrdY)
        result = interpol(xv1, velfga, xv2, velfgb, CrdX)
        Vel_V = result * (1.0 + 0.36 * veglngt / flngt) / (np.log((Altura + 0.36 * veglngt) / (0.13 * veglngt))) * (
                np.log((flngt / veglngt + 0.36) / 0.13) - 1.0)

        # ------ Calculation of the W wind speed -------
        velfga = interpol(yv1, W[iv, jv, 1], yv2, W[iv, jv + 1, 1], CrdY)
        velfgb = interpol(yv1, W[iv + 1, jv, 1], yv2, W[iv + 1, jv + 1, 1], CrdY)
        result = interpol(xv1, velfga, xv2, velfgb, CrdX)

        Vel_W = result * (1.0 + 0.36 * veglngt / flngt) / (np.log((Altura + 0.36 * veglngt) / (0.13 * veglngt))) * (
                np.log((flngt / veglngt + 0.36) / 0.13) - 1.0)

        return Vel_U, Vel_V, Vel_W

    def interpol(self, x1, y1, x2, y2, x):
        return y1 + (y2 - y1) / (0.0000000001 + x2 - x1) * (x - x1)
