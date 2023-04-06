from typing import TextIO

import numpy as np

from config.business_processor import BusinessConfigProcessor
from config.environment.environment import WindIndex, Environment, NodesNumber, WindMeshSpacing


class LoadWind(BusinessConfigProcessor):

    def process(self, file: TextIO):

        wind_model_name = file.readline().strip()

        if wind_model_name[:6] == "Canyon":
            # Canyon
            nodes = list(map(float, file.readline().split()))
            NodesNumber.NiWnd, NodesNumber.NjWnd, NodesNumber.NkWnd = list(map(int, nodes[:3]))
            NodesNumber.AuxSingle = nodes[3]
            self.wind_matrix_dimensioning()
            WindIndex.X_Wnd[0], WindIndex.Y_Wnd[0], WindMeshSpacing.DxWnd = map(float, file.readline().split(","))
            WindMeshSpacing.DyWnd = WindMeshSpacing.DxWnd
        else:
            # Nuatmos
            nodes = map(int, file.readline().split())
            NodesNumber.NiWnd, NodesNumber.NjWnd, NodesNumber.NkWnd = nodes
            self.wind_matrix_dimensioning()

        # The index 0 corresponds to the fictitious nodes.
        # Thus, the mesh defined initially corresponds, in i, from 1 to NiWnd
        NodesNumber.NiWnd -= 1
        NodesNumber.NjWnd -= 1
        NodesNumber.NkWnd -= 1

        if wind_model_name[:6] == "Canyon":
            # Canyon
            for i in range(NodesNumber.NiWnd):
                for j in range(NodesNumber.NjWnd):
                    for k in range(NodesNumber.NkWnd):
                        WindIndex.U[i][j][k], WindIndex.V[i][j][k], WindIndex.W[i][j][k], teijk, edijk, tijk, pijk, \
                            WindIndex.Z_Wnd[i][j][k] = map(float, file.readline().split(","))

            for i in range(1, NodesNumber.NiWnd):
                WindIndex.X_Wnd[i] = WindIndex.X_Wnd[0] + i * WindMeshSpacing.DxWnd

            for j in range(1, NodesNumber.NjWnd):
                WindIndex.Y_Wnd[j] = WindIndex.Y_Wnd[0] + j * WindMeshSpacing.DxWnd

        else:
            # -- Nuatmos
            for i in range(NodesNumber.NiWnd):
                for j in range(NodesNumber.NjWnd):
                    for k in range(NodesNumber.NkWnd):
                        WindIndex.X_Wnd[i], WindIndex.Y_Wnd[j], WindIndex.Z_Wnd[i][j][k] = \
                            map(float, file.readline().split())
                        WindIndex.U[i][j][k], WindIndex.V[i][j][k], WindIndex.W[i][j][k] = \
                            map(float, file.readline().split())

            WindMeshSpacing.DxWnd = WindIndex.X_Wnd[2] - WindIndex.X_Wnd[1]
            WindMeshSpacing.DyWnd = WindIndex.Y_Wnd[2] - WindIndex.Y_Wnd[1]

        return self.final_calculations()

    @staticmethod
    def wind_matrix_dimensioning():
        # Create and initialize arrays
        WindIndex.X_Wnd = np.zeros(NodesNumber.NiWnd + 1)
        WindIndex.Y_Wnd = np.zeros(NodesNumber.NjWnd + 1)
        WindIndex.Z_Wnd = np.zeros((NodesNumber.NiWnd + 1, NodesNumber.NjWnd + 1, NodesNumber.NkWnd + 1))

        WindIndex.XX_Wnd = np.zeros(NodesNumber.NiWnd)
        WindIndex.YY_Wnd = np.zeros(NodesNumber.NjWnd)
        WindIndex.ZZ_Wnd = np.zeros((NodesNumber.NiWnd + 1, NodesNumber.NjWnd + 1, NodesNumber.NkWnd + 1))

        WindIndex.U = np.zeros((NodesNumber.NiWnd + 1, NodesNumber.NjWnd + 1, NodesNumber.NkWnd + 1))
        WindIndex.Y = np.zeros((NodesNumber.NiWnd + 1, NodesNumber.NjWnd + 1, NodesNumber.NkWnd + 1))
        WindIndex.W = np.zeros((NodesNumber.NiWnd + 1, NodesNumber.NjWnd + 1, NodesNumber.NkWnd + 1))

    @staticmethod
    def final_calculations():
        for i in range(NodesNumber.NiWnd - 1):
            WindIndex.XX_Wnd[i] = 0.5 * (WindIndex.X_Wnd[i] + WindIndex.X_Wnd[i + 1])

        for j in range(NodesNumber.NjWnd - 1):
            WindIndex.YY_Wnd[j] = 0.5 * (WindIndex.Y_Wnd[j] + WindIndex.Y_Wnd[j + 1])

        Environment.X_IniWnd = WindIndex.X_Wnd[1]
        Environment.Y_IniWnd = WindIndex.Y_Wnd[1]
        Environment.X_FinWnd = WindIndex.X_Wnd[NodesNumber.NiWnd]
        Environment.Y_FinWnd = WindIndex.Y_Wnd[NodesNumber.NjWnd]

        for i in range(NodesNumber.NiWnd - 1):
            for j in range(NodesNumber.NjWnd - 1):
                for k in range(NodesNumber.NkWnd - 1):
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

    @staticmethod
    def interpol_wind_speed():
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

    @staticmethod
    def interpol(x1, y1, x2, y2, x):
        return y1 + (y2 - y1) / (0.0000000001 + x2 - x1) * (x - x1)


