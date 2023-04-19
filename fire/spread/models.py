import math

import numpy as np

from fire.element import CombustibleElement


def get_rate_of_spread(element1: CombustibleElement, element2: CombustibleElement):
    """
    TODO: For now, same direction in both elements, needs to be checked later and updated. Approved by professor Lopes
    :return: Return harmonic mean for element's rate of spread in direction gamma (ft/min)
    """
    angle_between_elements = get_angle(element1, element2)

    rate1 = get_rate_of_spread_direction_gamma(element1, angle_between_elements)
    rate2 = get_rate_of_spread_direction_gamma(element2, angle_between_elements)

    return (2 * rate1 * rate2) / (rate1 + rate2)


def get_rate_of_spread_direction_gamma(element: CombustibleElement, angle_between_elements: float) -> float:
    """
    Rate of spread in direction γ (gamma) from the ignition point relative to max spread direction (radian)
    References: Rγ in the equation of Patricia etc. 2018, p:86
    Max_spread_direction is relative to upslope. We're summing it with aspect to find the absolute_max_spread_direction.
    :return: Return rate of spread in direction gamma
    """
    absolute_max_spread_direction = element.aspect + element.max_spread_direction
    absolute_angle_between_elements = angle_between_elements - absolute_max_spread_direction

    rate_of_spread_gamma = element.rate_of_spread_heading_fire * (1 - element.e) / (1 - element.e * math.cos(absolute_angle_between_elements))
    return rate_of_spread_gamma


def get_angle(element1: CombustibleElement, element2: CombustibleElement) -> float:
    """
    :return: Return the angle between the elements (measured in radians)
    """
    latitude_distance = element2.location.latitude - element1.location.latitude
    longitude_distance = element2.location.longitude - element1.location.longitude

    return math.atan2(latitude_distance, longitude_distance)


def get_fireline_intensity(element1: CombustibleElement, element2: CombustibleElement) -> float:
    """
    :return: Return fireline intensity for directions other than head fire (Btu/ft/s)
    """
    return element1.heat_per_unit_area * get_rate_of_spread_direction_psi(element1, element2) / 60


def get_rate_of_spread_direction_psi(element1: CombustibleElement, element2: CombustibleElement) -> float:
    """
    Rate of spread in direction of normal to fire perimeter ψ (psi) at the point associated with the direction
    from ignition point γ (gamma), directions relative to direction of maximum spread (ft/min)
    TODO: psi direction is relative to max spread direction! Consider it!
    :return:
    """

    angle_between_elements = get_angle(element1, element2)
    absolute_max_spread_direction = element1.aspect + element1.max_spread_direction

    # Direction from the ignition point (y) relative to max spread direction (radian)
    absolute_neighbor_direction = angle_between_elements - absolute_max_spread_direction

    theta = np.arccos((element1.h * math.cos(absolute_neighbor_direction) *
                       (element1.h ** 2 * math.cos(absolute_neighbor_direction) ** 2 +
                        (element1.f ** 2 - element1.g ** 2) * math.sin(absolute_neighbor_direction) ** 2) ** 0.5 -
                       element1.f * element1.g * math.sin(absolute_neighbor_direction) ** 2) /
                      (element1.h ** 2 * math.cos(absolute_neighbor_direction) ** 2 + element1.f ** 2 *
                       math.sin(absolute_neighbor_direction) ** 2))

    r_psi = (element1.rate_of_spread_heading_fire * element1.h * (element1.g * math.cos(theta) + element1.f)) / \
            (element1.h ** 2 * math.cos(theta) ** 2 + element1.f ** 2 * math.sin(theta) ** 2) ** 0.5
    return r_psi
