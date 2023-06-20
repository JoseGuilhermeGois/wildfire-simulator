from dataclasses import dataclass


@dataclass
class Elevation:
    up_slope_direction_list: list[list[float]]
    slope_steepness_list: list[list[float]]
