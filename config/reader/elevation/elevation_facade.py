import math

import numpy as np

from config import Landscape
from fire import BaseTerrainTopographyFacade


class ElevationFacade:

    @staticmethod
    def calc_slope(number_of_columns, number_of_rows, element_size, elevation_distribution):
        up_slope_direction_list = np.zeros((number_of_rows + 1, number_of_columns + 1))
        slope_steepness_list = np.zeros((number_of_rows + 1, number_of_columns + 1))

        for j in range(2, number_of_rows - 1):
            for i in range(2, number_of_columns - 1):
                xxi = 2 * element_size
                yet = 2 * element_size
                zxi = float(elevation_distribution[j][i + 1]) - float(elevation_distribution[j][i - 1])
                zet = float(elevation_distribution[j + 1][i]) - float(elevation_distribution[j - 1][i])
                zex = -zxi / xxi
                zey = -zet / yet
                zez = 1

                # DirMaxDec() - Direction of the slope, measured from xx axis
                up_slope_direction_list[j][i] = 10 * ((BaseTerrainTopographyFacade.calc_angle(zey, zex) - math.pi) *
                                                      180 / math.pi)

                # AngIncl() - Slope steepness angle (always positive)
                zex2 = zex * zex
                zey2 = zey * zey
                zez2 = 1

                if zex2 == 0 and zey2 == 0:
                    slope_steepness_list[j][i] = 0
                else:
                    numerator = (zex2 + zey2) / (math.sqrt(zex2 + zey2 + zez2) * math.sqrt(zex2 + zey2))
                    slope_steepness_list[j][i] = 10 * (0.5 * math.pi - math.acos(numerator)) * 180 / math.pi

        # Corners for slope and dmc
        for i in range(2, number_of_columns - 1):
            up_slope_direction_list[1][i] = up_slope_direction_list[2][i]
            up_slope_direction_list[number_of_rows][i] = up_slope_direction_list[number_of_rows - 1][i]

            slope_steepness_list[1][i] = slope_steepness_list[2][i]
            slope_steepness_list[number_of_rows][i] = slope_steepness_list[number_of_rows - 1][i]

        for j in range(2, number_of_rows - 1):
            up_slope_direction_list[j][1] = up_slope_direction_list[j][2]
            up_slope_direction_list[j][number_of_columns] = up_slope_direction_list[j][number_of_columns - 1]

            slope_steepness_list[j][1] = slope_steepness_list[j][2]
            slope_steepness_list[j][number_of_columns] = slope_steepness_list[j][number_of_columns - 1]

        up_slope_direction_list[1][1] = up_slope_direction_list[2][2]
        up_slope_direction_list[number_of_rows][1] = up_slope_direction_list[number_of_rows - 1][2]
        up_slope_direction_list[1][number_of_columns] = up_slope_direction_list[2][number_of_columns - 1]
        up_slope_direction_list[number_of_rows][number_of_columns] = up_slope_direction_list[
            number_of_rows - 1][number_of_columns - 1]

        slope_steepness_list[1][1] = slope_steepness_list[2][2]
        slope_steepness_list[number_of_rows][1] = slope_steepness_list[number_of_rows - 1][2]
        slope_steepness_list[1][number_of_columns] = slope_steepness_list[2][number_of_columns - 1]
        slope_steepness_list[number_of_rows][number_of_columns] = slope_steepness_list[number_of_rows - 1][
            number_of_columns - 1]

        return up_slope_direction_list, slope_steepness_list
