import math
import sys
from typing import Protocol

import numpy as np

from config.environment import Environment
from config.fuel import Fuel
from fire.element import CombustibleElement, IncombustibleElement, Element, State
from config.landscape import Landscape, Location
from config.elevation import Elevation


class TerrainTopographyFacade(Protocol):

    def get_element(self, longitude: int, latitude: int) -> Element:
        ...


class BaseTerrainTopographyFacade(TerrainTopographyFacade):

    def __init__(self, landscape: Landscape, fuel_models: dict[str, Fuel], environment: list[list[Environment]],
                 elevation: Elevation):
        self.landscape: Landscape = landscape
        self.environment: list[list[Environment]] = environment
        self.fuel_models: dict[str, Fuel] = fuel_models
        self.elevation: Elevation = elevation
        self.element_size: int = landscape.element_size

    def get_element(self, longitude: int, latitude: int) -> Element:
        fuel_model_id: str = self.landscape.fuel_model_distribution[longitude][latitude]
        fuel: Fuel = self.fuel_models[fuel_model_id]
        environment_values: Environment = self.environment[longitude][latitude]
        elevation_horizontal_neighbours = (float(self.elevation.elevation_distribution[longitude - 1][latitude]),
                                           float(self.elevation.elevation_distribution[longitude + 1][latitude])
                                           )
        elevation_vertical_neighbours = (float(self.elevation.elevation_distribution[longitude][latitude - 1]),
                                         float(self.elevation.elevation_distribution[longitude][latitude + 1])
                                         )

        if fuel is None or fuel_model_id == '0':
            return IncombustibleElement()

        return CombustibleElement(
            location=Location(latitude, longitude),
            spread_time=sys.maxsize,
            state=State.FLAMMABLE,
            slope_steepness=self.calculate_slope_steepness(elevation_horizontal_neighbours,
                                                           elevation_vertical_neighbours),
            upslope_direction=self.calculate_upslope_direction(elevation_horizontal_neighbours,
                                                               elevation_vertical_neighbours),
            r_0=self.rate_of_spread_without_wind_slope(fuel),
            r_wind_up_slope=self.rate_of_spread_wind_upslope(fuel, elevation_horizontal_neighbours,
                                                             elevation_vertical_neighbours, environment_values),
            residence_time=self.residence_time(fuel),
            heat_per_unit_area=self.heat_per_unit_area(fuel),
            fireline_intensity_heading_fire=self.fireline_intensity_heading_fire(fuel, elevation_horizontal_neighbours,
                                                                                 elevation_vertical_neighbours,
                                                                                 environment_values),
            flame_length=self.flame_length(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours,
                                           environment_values),
            rate_of_spread_heading_fire=self.rate_of_spread_heading_fire(fuel, elevation_horizontal_neighbours,
                                                                            elevation_vertical_neighbours,
                                                                            environment_values),
            max_spread_direction=self.max_spread_direction(fuel, elevation_horizontal_neighbours,
                                                           elevation_vertical_neighbours, environment_values),
            effective_wind_speed=self.effective_wind_speed(fuel, elevation_horizontal_neighbours,
                                                           elevation_vertical_neighbours, environment_values),
            rate_of_spread_backing_fire=self.rate_of_spread_backing_fire(fuel, elevation_horizontal_neighbours,
                                                                         elevation_vertical_neighbours,
                                                                         environment_values),
            e=self.eccentricity(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours,
                                environment_values),
            f=self.ellipse_dimensions(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours,
                                      environment_values)[0],
            g=self.ellipse_dimensions(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours,
                                      environment_values)[1],
            h=self.ellipse_dimensions(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours,
                                      environment_values)[2]
        )

    # Unit: fraction
    def zex(self, elevation_horizontal_neighbours) -> float:
        """horizontal value for slope steepness and aspect calculation"""
        return -(elevation_horizontal_neighbours[1] - elevation_horizontal_neighbours[0]) / 2 * (self.element_size /
                                                                                                 3.281)

    # Unit: fraction
    def zey(self, elevation_vertical_neighbours) -> float:
        """vertical value for slope steepness and aspect calculation"""
        return -(elevation_vertical_neighbours[1] - elevation_vertical_neighbours[0]) / 2 * (self.element_size /
                                                                                               3.281)

    # Unit: radian (fraction)
    def calculate_slope_steepness(self, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        """Calculation of slope steepness from elevation"""
        zex2 = self.zex(elevation_horizontal_neighbours) ** 2
        zey2 = self.zey(elevation_vertical_neighbours) ** 2
        if zex2 == 0 and zey2 == 0:  # TODO: check this clause later
            return 0
        else:
            return math.sqrt(zex2 + zey2)

    # Unit: degree (azimuth)
    def calculate_aspect(self, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        """Calculation of aspect from elevation"""
        if self.zex(elevation_horizontal_neighbours) == 0 or self.zey(elevation_vertical_neighbours) == 0:  # TODO: check this clause later
            aspect = 0
        else:
            aspect = 180 - np.arctan(self.zey(elevation_vertical_neighbours) / self.zex(elevation_horizontal_neighbours) + 90 * (self.zex(elevation_horizontal_neighbours) / abs(self.zex(elevation_horizontal_neighbours))))

        aspect = 180 - aspect  # to fix coordination difference between numpy array and cardinal direction
        aspect += 360 if aspect < 0 else 0  # to get positive degree
        return aspect

    # Unit: degree
    def calculate_upslope_direction(self, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        """Up-slope direction is the opposite of aspect"""
        upslope_direction = self.calculate_aspect(elevation_horizontal_neighbours, elevation_vertical_neighbours) - 180
        upslope_direction += 360 if upslope_direction < 0 else 0
        return upslope_direction

    # Unit: m^2
    # area_ij: [[dead_1h, dead_10h, dead_100h], [live_herb, live_woody]]
    def area_ij(self, fuel) -> list[list[float]]:
        """Mean total surface area per unit fuel cell of each size class within each category"""
        area_ij = [[0] * 3, [0] * 2]
        for i in range(len(area_ij)):
            for j in range(len(area_ij[i])):
                area_ij[i][j] = fuel.sav_ratio[i][j] * fuel.load[i][j] / fuel.particle_density
        return area_ij

    # Unit: m^2
    # area_i: (dead, live)
    def area_i(self, fuel) -> tuple[float, float]:
        """Total mean surface area of the live and dead categories separately"""
        return sum(self.area_ij(fuel)[0]), sum(self.area_ij(fuel)[1])

    # Unit: m^2
    def area_t(self, fuel) -> float:
        """Total mean surface area of the live and dead categories all together"""
        return sum(self.area_i(fuel))

    # Unit: fraction
    # factor_ij: [[dead_1h, dead_10h, dead_100h], [live_herb, live_woody]]
    def factor_ij(self, fuel) -> list[list[float]]:
        """Weighting factor for characteristic dead and live heat content, effective mineral content,
        moisture content and surface-area-to-volume ratio"""
        factor_ij = [[0] * 3, [0] * 2]
        for i in range(len(factor_ij)):
            for j in range(len(factor_ij[i])):
                factor_ij[i][j] = self.area_ij(fuel)[i][j] / self.area_i(fuel)[i] \
                    if self.area_i(fuel)[i] > 0 else 0
        return factor_ij

    # Unit: fraction
    # factor_i: (dead, live)
    def factor_i(self, fuel) -> tuple[float, float]:
        """Weighting factor for characteristic fuel bed surface-area-to-volume ratio"""
        return (self.area_i(fuel)[0] / self.area_t(fuel),
                    self.area_i(fuel)[1] / self.area_t(fuel))

    # Unit: lb/ft^2
    # net_load: [[dead_1h, dead_10h, dead_100h], [live_herb, live_woody]]
    def net_load(self, fuel) -> list[list[float]]:
        """Net load of each size class within each category"""
        net_load = [[0] * 3, [0] * 2]
        for i in range(len(net_load)):
            for j in range(len(net_load[i])):
                net_load[i][j] = fuel.load[i][j] * (1 - fuel.total_mineral_content)
        return net_load

    # Unit: lb/ft^2
    # omg0d in Firestation
    def net_load_dead(self, fuel) -> float:
        """Net load of dead fuel category"""
        return self.factor_ij(fuel)[0][0] * self.net_load(fuel)[0][0] + \
                        self.factor_ij(fuel)[0][1] * self.net_load(fuel)[0][1] + \
                        self.factor_ij(fuel)[0][2] * self.net_load(fuel)[0][2]

    # Unit: lb/ft^2
    def net_load_live(self, fuel) -> float:
        """Net loaf of live fuel category"""
        return self.net_load(fuel)[1][0] + self.net_load(fuel)[1][1]

    # Unit: Btu/lb
    def heat_content_i(self, fuel) -> tuple[float, float]:
        return (sum(self.factor_ij(fuel)[0]) * fuel.heat_content[0],
                sum(self.factor_ij(fuel)[1]) * fuel.heat_content[1])

    # Unit: fraction
    def effective_mineral_content_i(self, fuel) -> tuple[float, float]:
        if sum(self.factor_ij(fuel)[1]) > 0:
            return (sum(self.factor_ij(fuel)[0]) * fuel.effective_mineral_content,
                                           sum(self.factor_ij(fuel)[1]) * fuel.effective_mineral_content)
        else:
            return sum(self.factor_ij(fuel)[0]) * fuel.effective_mineral_content, 0

    # Unit: fraction
    def mineral_damping_coff_i(self, fuel) -> tuple[float, float]:
        if sum(self.factor_ij(fuel)[1]) > 0:
            return (min(0.174 * self.effective_mineral_content_i(fuel)[0] ** (-0.19), 1),
                                      min(0.174 * self.effective_mineral_content_i(fuel)[1] ** (-0.19), 1))
        else:
            return min(0.174 * self.effective_mineral_content_i(fuel)[0] ** (-0.19), 1), 0

    # Unit: fraction
    # aamfd in Firestation
    def moisture_content_i(self, fuel) -> tuple[float, float]:

        return (self.factor_ij(fuel)[0][0] * fuel.moisture_content[0][0] +
                              self.factor_ij(fuel)[0][1] * fuel.moisture_content[0][1] +
                              self.factor_ij(fuel)[0][2] * fuel.moisture_content[0][2],
                              self.factor_ij(fuel)[1][0] * fuel.moisture_content[1][0] +
                              self.factor_ij(fuel)[1][1] * fuel.moisture_content[1][1])

    # Unit: fraction
    def dead_to_live_load_ratio(self, fuel) -> float:
        dead_to_live_load_ratio_num = 0
        dead_to_live_load_ratio_den = 0
        for j in range(len(fuel.sav_ratio[0])):
            if fuel.sav_ratio[0][j] > 0:
                dead_to_live_load_ratio_num += fuel.load[0][j] * np.exp(-138 / (fuel.sav_ratio[0][j]))

        for j in range(len(fuel.sav_ratio[1])):
            if fuel.sav_ratio[1][j] > 0:
                dead_to_live_load_ratio_den += fuel.load[1][j] * np.exp(-500 / (fuel.sav_ratio[1][j]))

        if dead_to_live_load_ratio_den > 0:
            return dead_to_live_load_ratio_num / dead_to_live_load_ratio_den
        else:
            return 0

    # Unit: fraction
    def fine_dead_fuel_moisture(self, fuel) -> float:
        fine_dead_fuel_moisture_num = 0
        fine_dead_fuel_moisture_den = 0
        for j in range(len(fuel.sav_ratio[0])):
            if fuel.sav_ratio[0][j] > 0:
                fine_dead_fuel_moisture_num += fuel.moisture_content[0][j] * fuel.load[0][j] * \
                                               np.exp(-138 / (fuel.sav_ratio[0][j]))
        for j in range(len(fuel.sav_ratio[0])):
            if fuel.sav_ratio[0][j] > 0:
                fine_dead_fuel_moisture_den += fuel.load[0][j] * np.exp(-138 / (fuel.sav_ratio[0][j]))

        return fine_dead_fuel_moisture_num / fine_dead_fuel_moisture_den if fine_dead_fuel_moisture_den > 0 else 0

    # Unit: fraction
    def moisture_extinction_live(self, fuel) -> float:
        return max(2.9 * self.dead_to_live_load_ratio(fuel) * (1 - (self.fine_dead_fuel_moisture(fuel) /
                                    fuel.moisture_extinction_dead)) - 0.226, fuel.moisture_extinction_dead)

    # Unit: fraction
    def tau_m(self, fuel) -> list[float, float]:
        return [min(self.moisture_content_i(fuel)[0] / fuel.moisture_extinction_dead, 1),
                 min(self.moisture_content_i(fuel)[1] / self.moisture_extinction_live(fuel), 1)]

    # Unit: fraction
    def moisture_damping_coff_i(self, fuel) -> list[float, float]:
        moisture_damping_coff_i = [0] * 2
        for i in range(0, 2):
            moisture_damping_coff_i[i] = 1 - 2.59 * self.tau_m(fuel)[i] + 5.11 * (self.tau_m(fuel)[i] ** 2) - 3.52 * \
                                         (self.tau_m(fuel)[i] ** 3)
        return moisture_damping_coff_i

    # Unit: ft^-1
    # Sigmad(n) in Firestation
    def sav_ratio_i(self, fuel) -> tuple[float,float]:
        return (self.factor_ij(fuel)[0][0] * fuel.sav_ratio[0][0] + self.factor_ij(fuel)[0][1] * fuel.sav_ratio[0][1] +
                       self.factor_ij(fuel)[0][2] * fuel.sav_ratio[0][2], self.factor_ij(fuel)[1][0] * fuel.sav_ratio[1][0] +
                       self.factor_ij(fuel)[1][1] * fuel.sav_ratio[1][1])

    # Unit: ft^-1
    def sav_ratio_characteristic(self, fuel) -> float:
        return self.factor_i(fuel)[0] * self.sav_ratio_i(fuel)[0] + self.factor_i(fuel)[1] * self.sav_ratio_i(fuel)[1]

    # Unit: lb/ft^3
    def mean_bulk_density(self, fuel) -> float:
        return (sum(fuel.load[0]) + sum(fuel.load[1])) / fuel.bed_depth

    # Unit: fraction
    def mean_packing_ratio(self, fuel) -> float:
        return self.mean_bulk_density(fuel) / fuel.particle_density

    # Unit: fraction
    def optimum_packing_ratio(self, fuel) -> float:
        return 3.348 * self.sav_ratio_characteristic(fuel) ** (-0.8189)

    # Unit: fraction
    def relative_packing_ratio(self, fuel) -> float:
        return self.mean_packing_ratio(fuel) / self.optimum_packing_ratio(fuel)

    # Unit: fraction
    def propagating_flux_ratio(self, fuel) -> float:
        """Heat Source"""
        return (192 + 0.2595 * self.sav_ratio_characteristic(fuel)) ** (-1) * np.exp(
            (0.792 + 0.681 * self.sav_ratio_characteristic(fuel) ** 0.5) * (self.mean_packing_ratio(fuel) + 0.1))

    # Unit: min^-1
    def aaa(self, fuel) -> float:
        """A exponent"""
        return 133 * self.sav_ratio_characteristic(fuel) ** (-0.7913)

    # Unit: min^-1
    def maximum_reaction_velocity(self, fuel) -> float:
        return self.sav_ratio_characteristic(fuel) ** 1.5 * (495 + 0.0594 * self.sav_ratio_characteristic(fuel) ** 1.5) ** (-1)

    # Unit: min^-1
    def optimum_reaction_velocity(self, fuel) -> float:
        return self.maximum_reaction_velocity(fuel) * (self.relative_packing_ratio(fuel) ** self.aaa(fuel)) * \
                                    np.exp(self.aaa(fuel) * (1 - self.relative_packing_ratio(fuel)))

    # Unit: Btu/ft -min
    def reaction_intensity_dead(self, fuel) -> float:
        return self.net_load_dead(fuel) * self.heat_content_i(fuel)[0] * self.moisture_damping_coff_i(fuel)[0] * \
                self.mineral_damping_coff_i(fuel)[0]

    # Unit: Btu/ft -min
    def reaction_intensity_live(self, fuel) -> float:
        return self.net_load_live(fuel) * self.heat_content_i(fuel)[1] * self.moisture_damping_coff_i(fuel)[1] * \
                self.mineral_damping_coff_i(fuel)[1]

    # Unit: Btu/ft -min
    def reaction_intensity(self, fuel) -> float:
        return self.optimum_reaction_velocity(fuel) * (self.reaction_intensity_dead(fuel) + self.reaction_intensity_live(fuel))

    # Unit: Btu/lb
    def heat_of_pre_ignition_ij(self, fuel) -> list[list[float]]:
        heat_of_pre_ignition_ij = [[0] * 3, [0] * 2]
        for i in range(2):
            for j in range(len(heat_of_pre_ignition_ij[i])):
                heat_of_pre_ignition_ij[i][j] = 250 + 1116 * fuel.moisture_content[i][j]
        return heat_of_pre_ignition_ij

    # Unit: fraction
    def heat_number_tmp(self, fuel) -> list[float]:
        heat_number_tmp = [0] * 2
        for i in range(len(heat_number_tmp)):
            heat_number_tmp[i] = self.factor_i(fuel)[i] * sum([a * (np.exp(-138 / b) if b > 0 else 0) * c
                                                     for a, b, c in zip(self.factor_ij(fuel)[i],
                                                                        fuel.sav_ratio[i],
                                                                        self.heat_of_pre_ignition_ij(fuel)[i])])
        return heat_number_tmp

    # Unit: Btu/ft^3
    def heat_sink(self, fuel) -> float:
        """The denominator (heat sink) of ROS calculation is the heat required for ignition by the potential fuel"""
        return self.mean_bulk_density(fuel) * sum(self.heat_number_tmp(fuel))

    def wind_limit(self, fuel) -> float:
        """It is now recommended that a wind limit NOT be imposed (Andrews et al. 2013)"""
        # return min(self.wind_speed, 96.8 * self.reaction_intensity() ** 1 / 3)
        pass

    def ccc(self, fuel):
        return 7.47 * np.exp(-0.133 * (self.sav_ratio_characteristic(fuel) ** 0.55))

    def bbb(self, fuel):
        return 0.02526 * (self.sav_ratio_characteristic(fuel) ** 0.54)

    def eee(self, fuel):
        return 0.715 * np.exp(-0.000359 * self.sav_ratio_characteristic(fuel))

    def wind_factor(self, fuel, enviroment_values) -> float:
        return self.ccc(fuel) * (enviroment_values.wind_speed ** self.bbb(fuel)) * (self.relative_packing_ratio(fuel) ** (self.eee(fuel) * -1))

    def slope_factor(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        return 5.275 * (self.mean_packing_ratio(fuel) ** (-0.3)) * (self.calculate_slope_steepness(elevation_horizontal_neighbours, elevation_vertical_neighbours) ** 2)

    # Unit: Btu/ft^2-min
    def heat_source(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, environment_values) -> float:
        """The numerator (heat source) of ROS calculation is reaction intensity times propagating flux ratio
        increased by the influence of wind and slope"""
        return self.reaction_intensity(fuel) * self.propagating_flux_ratio(fuel) * (1 + self.wind_factor(fuel, environment_values) + self.slope_factor(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours))

    # Unit: ft/min
    def rate_of_spread_without_wind_slope(self, fuel) -> float:
        """Rate of spread with no wind and slope"""
        return self.reaction_intensity(fuel) * self.propagating_flux_ratio(fuel) / self.heat_sink(fuel)

    # Unit: ft/min
    def rate_of_spread_wind_upslope(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """Rate of spread if wind is up slope"""
        return self.heat_source(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) / self.heat_sink(fuel)

    # Unit: degree
    def wind_direction(self, enviroment_values) -> float:
        """The wind direction used in this model is the opposite of the given wind direction.
        Because the given wind direction is defined as the direction the wind is coming from.
        For this model, we need the direction the wind goes to."""
        wind_direction = enviroment_values.wind_direction - 180
        wind_direction += 360 if wind_direction < 0 else 0
        return wind_direction

    # Unit: degree
    def omega(self, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """Wind direction ω (omega) relative to up-slope"""
        omega = self.wind_direction(enviroment_values) - self.calculate_upslope_direction(elevation_horizontal_neighbours, elevation_vertical_neighbours)
        omega += 360 if omega < 0 else 0
        return omega

    def magnitude_of_slope_vector(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        """The slope vector has magnitude Ds and direction 0, Andrews 2018, page. 86"""
        return self.rate_of_spread_without_wind_slope(fuel) * self.slope_factor(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours)

    def magnitude_of_wind_vector(self, fuel, enviroment_values) -> float:
        """The wind vector has magnitude Dw in direction ω (omega) from up-slope, Andrews 2018, page. 86"""
        return self.rate_of_spread_without_wind_slope(fuel) * self.wind_factor(fuel, enviroment_values)

    def resultant_vector_x(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """The resultant vector is then (Ds + Dw cosω, Dw sinω), Andrews 2018, page. 85"""
        return self.magnitude_of_slope_vector(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours) + (self.magnitude_of_wind_vector(fuel, enviroment_values) *
                                                   math.cos(math.radians(self.omega(elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values))))

    def resultant_vector_y(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """The resultant vector is then (Ds + Dw cosω, Dw sinω), Andrews 2018, page. 85"""
        return self.magnitude_of_wind_vector(fuel, enviroment_values) * math.sin(math.radians(self.omega(elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values)))

    def magnitude_of_heading_fire_vector(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """Compound vector (Dh) of slope and wind vectors in direction α (alpha), Andrews 2018, page. 86"""
        return (self.resultant_vector_x(fuel, elevation_horizontal_neighbours,
                                        elevation_vertical_neighbours,
                                        enviroment_values) ** 2
                +

                self.resultant_vector_y(fuel,
                                        elevation_horizontal_neighbours,
                                        elevation_vertical_neighbours,
                                        enviroment_values) ** 2) ** 0.5

    # Unit: degree
    def max_spread_direction(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours,
                             enviroment_values) -> float:
        # Direction of max spread α (alpha) relative to up slope (radian)
        return math.asin(abs(self.resultant_vector_y(fuel,
                                                     elevation_horizontal_neighbours,
                                                     elevation_vertical_neighbours,
                                                     enviroment_values))
                         /

                         self.magnitude_of_heading_fire_vector(fuel,
                                                               elevation_horizontal_neighbours,
                                                               elevation_vertical_neighbours,
                                                               enviroment_values))

    # Unit: ft/min
    def rate_of_spread_heading_fire(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """Rate of spread of headinging fire in the direction of max spread"""
        return self.rate_of_spread_without_wind_slope(fuel) + self.magnitude_of_heading_fire_vector(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values)

    def effective_wind_factor(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """Effective wind factor in direction of maximum spread"""
        return (self.rate_of_spread_heading_fire(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) / self.rate_of_spread_without_wind_slope(fuel)) - 1

    # Unit: ft/min
    def effective_wind_speed(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """Effective wind speed in direction of maximum spread"""
        return ((self.effective_wind_factor(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) * (self.relative_packing_ratio(fuel) ** self.eee(fuel))) /
                self.ccc(fuel)) ** (1/self.bbb(fuel))

    def length_to_width_ratio(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values):
        """The ellipse length-to-width ratio (z) in the direction of max spread.
        Effective wind speed in mi/h in the formulation of Z."""
        # TODO: check this later
        return 1 + 0.25 * (self.effective_wind_speed(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) * 0.01136)

    def eccentricity(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values):
        """The eccentricity of the ellipse (e)"""
        return (((self.length_to_width_ratio(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) ** 2) - 1) ** 0.5) / self.length_to_width_ratio(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values)

    # Unit: ft/min
    def rate_of_spread_backing_fire(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values):
        return self.rate_of_spread_heading_fire(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) * (1 - self.eccentricity(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values)) / (1 + self.eccentricity(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values))

    # Unit: ft
    def fire_length(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values):
        return self.rate_of_spread_heading_fire(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) + self.rate_of_spread_backing_fire(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values)

    # Unit: ft
    def max_fire_width(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values):
        return self.length_to_width_ratio(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) / self.fire_length(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values)

    # Unit: ft
    def flanking_spread_distance(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values):
        return self.max_fire_width(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) / 2

    # Unit: ft
    def ellipse_dimensions(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values):
        fff = self.fire_length(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) / 2
        ggg = self.rate_of_spread_heading_fire(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) - fff
        hhh = self.flanking_spread_distance(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values)

        return fff, ggg, hhh

    # Unit: min
    def residence_time(self, fuel) -> float:
        return (384 / self.sav_ratio_characteristic(fuel)) * self.element_size ** 2

    # Unit: Btu/ft^2
    def heat_per_unit_area(self, fuel) -> float:
        return self.reaction_intensity(fuel) * self.residence_time(fuel)

    # Unit: Btu/ft/s
    def fireline_intensity_heading_fire(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        """Bayram's fire-line intensity. Fire-line intensity for the direction of heading fire"""
        # TODO: Implement this for the directions other than heading fire (check page 91),
        #  it is important because of crown fire implementation.
        return self.heat_per_unit_area(fuel) * self.rate_of_spread_heading_fire(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) / 60

    # Bayram's fire-line intensity for the directions other than heading fire (Btu/ft/s)
    #   TODO: Fill here!!

    # Unit: ft
    def flame_length(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        return 0.45 * self.fireline_intensity_heading_fire(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) ** 0.46

    # Unit: ft
    def flame_depth(self, fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) -> float:
        return self.rate_of_spread_heading_fire(fuel, elevation_horizontal_neighbours, elevation_vertical_neighbours, enviroment_values) * self.residence_time(fuel)


