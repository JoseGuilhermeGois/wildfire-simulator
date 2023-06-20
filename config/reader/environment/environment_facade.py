import math
from typing import TextIO, Protocol

import numpy as np

from config.reader.environment.environment import WindIndex


class LoadModels:

    @staticmethod
    def wind_matrix_dimensioning(wind_index):
        # Create and initialize arrays
        wind_index.X_Wnd = np.zeros(wind_index.NiWnd + 1)
        wind_index.Y_Wnd = np.zeros(wind_index.NjWnd + 1)
        wind_index.Z_Wnd = np.zeros((wind_index.NiWnd + 1, wind_index.NjWnd + 1,
                                    wind_index.NkWnd + 1))

        wind_index.XX_Wnd = np.zeros(wind_index.NiWnd)
        wind_index.YY_Wnd = np.zeros(wind_index.NjWnd)
        wind_index.ZZ_Wnd = np.zeros((wind_index.NiWnd + 1, wind_index.NjWnd + 1,
                                     wind_index.NkWnd + 1))

        wind_index.U = np.zeros((wind_index.NiWnd + 1, wind_index.NjWnd + 1,
                                wind_index.NkWnd + 1))
        wind_index.V = np.zeros((wind_index.NiWnd + 1, wind_index.NjWnd + 1,
                                wind_index.NkWnd + 1))
        wind_index.W = np.zeros((wind_index.NiWnd + 1, wind_index.NjWnd + 1,
                                wind_index.NkWnd + 1))

    @staticmethod
    def remove_fictitious_nodes(wind_index):
        # The index 0 corresponds to the fictitious nodes.
        # Thus, the mesh defined initially corresponds, in i, from 1 to NiWnd
        wind_index.NiWnd -= 1
        wind_index.NjWnd -= 1
        wind_index.NkWnd -= 1

    @staticmethod
    def final_calculations(wind_index):
        for i in range(wind_index.NiWnd - 1):
            wind_index.XX_Wnd[i] = 0.5 * (wind_index.X_Wnd[i] + wind_index.X_Wnd[i + 1])

        for j in range(wind_index.NjWnd - 1):
            wind_index.YY_Wnd[j] = 0.5 * (wind_index.Y_Wnd[j] + wind_index.Y_Wnd[j + 1])

        wind_index.X_IniWnd = wind_index.X_Wnd[1]
        wind_index.Y_IniWnd = wind_index.Y_Wnd[1]
        wind_index.X_FinWnd = wind_index.X_Wnd[wind_index.NiWnd]
        wind_index.Y_FinWnd = wind_index.Y_Wnd[wind_index.NjWnd]

        for i in range(wind_index.NiWnd - 1):
            for j in range(wind_index.NjWnd - 1):
                for k in range(wind_index.NkWnd - 1):
                    wind_index.ZZ_Wnd[i][j][k] = 0.125 * (wind_index.Z_Wnd[i][j][k] + wind_index.Z_Wnd[i + 1][j][k] +
                                                          wind_index.Z_Wnd[i + 1][j + 1][k] + wind_index.Z_Wnd[i][j + 1][
                                                             k] +
                                                          wind_index.Z_Wnd[i][j][k + 1] + wind_index.Z_Wnd[i + 1][j][
                                                             k + 1] +
                                                          wind_index.Z_Wnd[i + 1][j + 1][k + 1] +
                                                          wind_index.Z_Wnd[i][j + 1][k + 1] - 2 *
                                                          wind_index.Z_Wnd[i][j][1] - 2 * wind_index.Z_Wnd[i + 1][j][
                                                             1] - 2 *
                                                          wind_index.Z_Wnd[i + 1][j + 1][k] - 2 *
                                                          wind_index.Z_Wnd[i][j + 1][1])
