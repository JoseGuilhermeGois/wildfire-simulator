import math
import sys
from typing import Protocol

import numpy as np

from config.environment import Environment
from config.fuel import Fuel
from fire.element import CombustibleElement, IncombustibleElement, Element, State
from config.landscape import Landscape, Location


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

    # Unit: m^2
    # area_ij: [[dead_1h, dead_10h, dead_100h], [live_herb, live_woody]]
    def area_ij(self) -> list[list[float]]:
        """Mean total surface area per unit fuel cell of each size class within each category"""
        area_ij = [[0] * 3, [0] * 2]
        for i in range(len(area_ij)):
            for j in range(len(area_ij[i])):
                area_ij[i][j] = self.fuel.sav_ratio[i][j] * self.fuel.load[i][j] / self.fuel.particle_density
        return area_ij

    # Unit: m^2
    # area_i: (dead, live)
    def area_i(self) -> tuple[float, float]:
        """Total mean surface area of the live and dead categories separately"""
        return sum(self.area_ij()[0]), sum(self.area_ij()[1])

    # Unit: m^2
    def area_t(self) -> float:
        """Total mean surface area of the live and dead categories all together"""
        return sum(self.area_i())

    # Unit: fraction
    # factor_ij: [[dead_1h, dead_10h, dead_100h], [live_herb, live_woody]]
    def factor_ij(self) -> list[list[float]]:
        """Weighting factor for characteristic dead and live heat content, effective mineral content,
        moisture content and surface-area-to-volume ratio"""
        factor_ij = [[0] * 3, [0] * 2]
        for i in range(len(factor_ij)):
            for j in range(len(factor_ij[i])):
                factor_ij[i][j] = self.area_ij()[i][j] / self.area_i()[i] \
                    if self.area_i()[i] > 0 else 0
        return factor_ij

    # Unit: fraction
    # factor_i: (dead, live)
    def factor_i(self) -> tuple[float, float]:
        """Weighting factor for characteristic fuel bed surface-area-to-volume ratio"""
        return (self.area_i()[0] / self.area_t(),
                    self.area_i()[1] / self.area_t())

    # Unit: lb/ft^2
    # net_load: [[dead_1h, dead_10h, dead_100h], [live_herb, live_woody]]
    def net_load(self) -> list[list[float]]:
        """Net load of each size class within each category"""
        net_load = [[0] * 3, [0] * 2]
        for i in range(len(net_load)):
            for j in range(len(net_load[i])):
                net_load[i][j] = self.fuel.load[i][j] * (1 - self.fuel.total_mineral_content)
        return net_load

    # Unit: lb/ft^2
    # omg0d in Firestation
    def net_load_dead(self) -> float:
        """Net load of dead fuel category"""
        return self.factor_ij()[0][0] * self.net_load()[0][0] + \
                        self.factor_ij()[0][1] * self.net_load()[0][1] + \
                        self.factor_ij()[0][2] * self.net_load()[0][2]

    # Unit: lb/ft^2
    def net_load_live(self) -> float:
        """Net loaf of live fuel category"""
        return self.net_load()[1][0] + self.net_load()[1][1]

    # Unit: Btu/lb
    def heat_content_i(self) -> tuple[float, float]:
        return (sum(self.factor_ij()[0]) * self.fuel.heat_content[0],
                sum(self.factor_ij()[1]) * self.fuel.heat_content[1])

    # Unit: fraction
    def effective_mineral_content_i(self) -> tuple[float, float]:
        if sum(self.factor_ij()[1]) > 0:
            return (sum(self.factor_ij()[0]) * self.fuel.effective_mineral_content,
                                           sum(self.factor_ij()[1]) * self.fuel.effective_mineral_content)
        else:
            return sum(self.factor_ij()[0]) * self.fuel.effective_mineral_content, 0

    # Unit: fraction
    def mineral_damping_coff_i(self) -> tuple[float, float]:
        if sum(self.factor_ij()[1]) > 0:
            return (min(0.174 * self.effective_mineral_content_i()[0] ** (-0.19), 1),
                                      min(0.174 * self.effective_mineral_content_i()[1] ** (-0.19), 1))
        else:
            return min(0.174 * self.effective_mineral_content_i()[0] ** (-0.19), 1), 0

    # Unit: fraction
    # aamfd in Firestation
    def moisture_content_i(self) -> tuple[float, float]:

        return (self.factor_ij()[0][0] * GridInfo.moisture_contents[0][0] +
                              self.factor_ij()[0][1] * GridInfo.moisture_contents[0][1] +
                              self.factor_ij()[0][2] * GridInfo.moisture_contents[0][2],
                              self.factor_ij()[1][0] * GridInfo.moisture_contents[1][0] +
                              self.factor_ij()[1][1] * GridInfo.moisture_contents[1][1])

    # Unit: fraction
    def dead_to_live_load_ratio(self) -> float:
        dead_to_live_load_ratio_num = 0
        dead_to_live_load_ratio_den = 0
        for j in range(len(self.fuel.sav_ratio[0])):
            if self.fuel.sav_ratio[0][j] > 0:
                dead_to_live_load_ratio_num += self.fuel.load[0][j] * np.exp(-138 / (self.fuel.sav_ratio[0][j]))

        for j in range(len(self.fuel.sav_ratio[1])):
            if self.fuel.sav_ratio[1][j] > 0:
                dead_to_live_load_ratio_den += self.fuel.load[1][j] * np.exp(-500 / (self.fuel.sav_ratio[1][j]))

        if dead_to_live_load_ratio_den > 0:
            return dead_to_live_load_ratio_num / dead_to_live_load_ratio_den
        else:
            return 0

    # Unit: fraction
    def fine_dead_fuel_moisture(self) -> float:
        fine_dead_fuel_moisture_num = 0
        fine_dead_fuel_moisture_den = 0
        for j in range(len(self.fuel.sav_ratio[0])):
            if self.fuel.sav_ratio[0][j] > 0:
                fine_dead_fuel_moisture_num += GridInfo.moisture_content[0][j] * self.fuel.load[0][j] * \
                                               np.exp(-138 / (self.fuel.sav_ratio[0][j]))
        for j in range(len(self.fuel.sav_ratio[0])):
            if self.fuel.sav_ratio[0][j] > 0:
                fine_dead_fuel_moisture_den += self.fuel.load[0][j] * np.exp(-138 / (self.fuel.sav_ratio[0][j]))

        return fine_dead_fuel_moisture_num / fine_dead_fuel_moisture_den if fine_dead_fuel_moisture_den > 0 else 0

    # Unit: fraction
    def moisture_extinction_live(self) -> float:
        return max(2.9 * self.dead_to_live_load_ratio() * (1 - (self.fine_dead_fuel_moisture() /
                                    self.fuel.moisture_extinction_dead)) - 0.226, self.fuel.moisture_extinction_dead)

    # Unit: fraction
    def tau_m(self) -> list[float,float]:
        return [min(GridInfo.moisture_content_i()[0] / self.fuel.moisture_extinction_dead, 1),
                 min(GridInfo.moisture_content_i()[1] / self.moisture_extinction_live(), 1)]

    # Unit: fraction
    def moisture_damping_coff_i(self) -> list[float, float]:
        moisture_damping_coff_i = [0] * 2
        for i in range(0, 2):
            moisture_damping_coff_i[i] = 1 - 2.59 * self.tau_m()[i] + 5.11 * (self.tau_m()[i] ** 2) - 3.52 * \
                                         (self.tau_m()[i] ** 3)
        return moisture_damping_coff_i

    # Unit: ft^-1
    # Sigmad(n) in Firestation
    def sav_ratio_i(self) -> tuple[float,float]:
        return (self.factor_ij()[0][0] * self.fuel.sav_ratio[0][0] + self.factor_ij()[0][1] * self.fuel.sav_ratio[0][1] +
                       self.factor_ij()[0][2] * self.fuel.sav_ratio[0][2], self.factor_ij()[1][0] * self.fuel.sav_ratio[1][0] +
                       self.factor_ij()[1][1] * self.fuel.sav_ratio[1][1])

    # Unit: ft^-1
    def sav_ratio_characteristic(self) -> float:
        return self.factor_i()[0] * self.sav_ratio_i()[0] + self.factor_i()[1] * self.sav_ratio_i()[1]

    # Unit: lb/ft^3
    def mean_bulk_density(self) -> float:
        return (sum(self.fuel.load[0]) + sum(self.fuel.load[1])) / self.fuel.bed_depth

    # Unit: fraction
    def mean_packing_ratio(self) -> float:
        return self.mean_bulk_density() / self.fuel.particle_density

    # Unit: fraction
    def optimum_packing_ratio(self) -> float:
        return 3.348 * self.sav_ratio_characteristic() ** (-0.8189)

    # Unit: fraction
    def relative_packing_ratio(self) -> float:
        return self.mean_packing_ratio() / self.optimum_packing_ratio()

    # Unit: fraction
    def propagating_flux_ratio(self) -> float:
        """Heat Source"""
        return (192 + 0.2595 * self.sav_ratio_characteristic()) ** (-1) * np.exp(
            (0.792 + 0.681 * self.sav_ratio_characteristic() ** 0.5) * (self.mean_packing_ratio() + 0.1))

    # Unit: min^-1
    def aaa(self) -> float:
        """A exponent"""
        return 133 * self.sav_ratio_characteristic() ** (-0.7913)

    # Unit: min^-1
    def maximum_reaction_velocity(self) -> float:
        return self.sav_ratio_characteristic() ** 1.5 * (495 + 0.0594 * self.sav_ratio_characteristic() ** 1.5) ** (-1)

    # Unit: min^-1
    def optimum_reaction_velocity(self) -> float:
        return self.maximum_reaction_velocity() * (self.relative_packing_ratio() ** self.aaa()) * \
                                    np.exp(self.aaa() * (1 - self.relative_packing_ratio()))

    # Unit: Btu/ft -min
    def reaction_intensity_dead(self) -> float:
        return self.net_load_dead() * self.heat_content_i()[0] * self.moisture_damping_coff_i()[0] * \
                self.mineral_damping_coff_i()[0]

    # Unit: Btu/ft -min
    def reaction_intensity_live(self) -> float:
        return self.net_load_live() * self.heat_content_i()[1] * self.moisture_damping_coff_i()[1] * \
                self.mineral_damping_coff_i()[1]

    # Unit: Btu/ft -min
    def reaction_intensity(self) -> float:
        return self.optimum_reaction_velocity() * (self.reaction_intensity_dead() + self.reaction_intensity_live())

    # Unit: Btu/lb
    def heat_of_pre_ignition_ij(self) -> list[list[float]]:
        heat_of_pre_ignition_ij = [[0] * 3, [0] * 2]
        for i in range(2):
            for j in range(len(heat_of_pre_ignition_ij[i])):
                heat_of_pre_ignition_ij[i][j] = 250 + 1116 * GridInfo.moisture_contents[i][j]
        return heat_of_pre_ignition_ij

    # Unit: fraction
    def heat_number_tmp(self) -> list[float]:
        heat_number_tmp = [0] * 2
        for i in range(len(heat_number_tmp)):
            heat_number_tmp[i] = self.factor_i()[i] * sum([a * (np.exp(-138 / b) if b > 0 else 0) * c
                                                     for a, b, c in zip(self.factor_ij()[i],
                                                                        self.fuel.sav_ratio[i],
                                                                        self.heat_of_pre_ignition_ij()[i])])
        return heat_number_tmp

    # Unit: Btu/ft^3
    def heat_sink(self) -> float:
        """The denominator (heat sink) of ROS calculation is the heat required for ignition by the potential fuel"""
        return self.mean_bulk_density() * sum(self.heat_number_tmp())

    def wind_limit(self) -> float:
        """It is now recommended that a wind limit NOT be imposed (Andrews et al. 2013)"""
        # return min(self.wind_speed, 96.8 * self.reaction_intensity() ** 1 / 3)
        pass

    def ccc(self):
        return 7.47 * np.exp(-0.133 * (self.sav_ratio_characteristic() ** 0.55))

    def bbb(self):
        return 0.02526 * (self.sav_ratio_characteristic() ** 0.54)

    def eee(self):
        return 0.715 * np.exp(-0.000359 * self.sav_ratio_characteristic())

    def wind_factor(self) -> float:
        return self.ccc() * (self.wind_speed ** self.bbb()) * (self.relative_packing_ratio() ** (self.eee() * -1))

    def slope_factor(self) -> float:
        return 5.275 * (self.mean_packing_ratio() ** (-0.3)) * (self.slope_steepness ** 2)

    # Unit: Btu/ft^2-min
    def heat_source(self) -> float:
        """The numerator (heat source) of ROS calculation is reaction intensity times propagating flux ratio
        increased by the influence of wind and slope"""
        return self.reaction_intensity() * self.propagating_flux_ratio() * (1 + self.wind_factor() + self.slope_factor())

    # Unit: ft/min
    def rate_of_spread_without_wind_slope(self) -> float:
        """Rate of spread with no wind and slope"""
        return self.reaction_intensity() * self.propagating_flux_ratio() / self.heat_sink()

    # Unit: ft/min
    def rate_of_spread_wind_upslope(self) -> float:
        """Rate of spread if wind is up slope"""
        return self.heat_source() / self.heat_sink()

    # Unit: degree
    def wind_direction(self) -> float:
        """The wind direction used in this model is the opposite of the given wind direction.
        Because the given wind direction is defined as the direction the wind is coming from.
        For this model, we need the direction the wind goes to."""
        wind_direction = self.wind_direction - 180
        wind_direction += 360 if wind_direction < 0 else 0
        return wind_direction

    # Unit: degree
    def omega(self) -> float:
        """Wind direction ω (omega) relative to up-slope"""
        omega = self.wind_direction() - self.upslope_direction
        omega += 360 if omega < 0 else 0
        return omega

    def magnitude_of_slope_vector(self) -> float:
        """The slope vector has magnitude Ds and direction 0, Andrews 2018, page. 86"""
        return self.rate_of_spread_without_wind_slope() * self.slope_factor()

    def magnitude_of_wind_vector(self) -> float:
        """The wind vector has magnitude Dw in direction ω (omega) from up-slope, Andrews 2018, page. 86"""
        return self.rate_of_spread_without_wind_slope() * self.wind_factor()

    def resultant_vector_x(self) -> float:
        """The resultant vector is then (Ds + Dw cosω, Dw sinω), Andrews 2018, page. 85"""
        return self.magnitude_of_slope_vector() + (self.magnitude_of_wind_vector() *
                                                   math.cos(math.radians(self.omega())))

    def resultant_vector_y(self) -> float:
        """The resultant vector is then (Ds + Dw cosω, Dw sinω), Andrews 2018, page. 85"""
        return self.magnitude_of_wind_vector() * math.sin(math.radians(self.omega()))

    def magnitude_of_heading_fire_vector(self) -> float:
        """Compound vector (Dh) of slope and wind vectors in direction α (alpha), Andrews 2018, page. 86"""
        return (self.resultant_vector_x() ** 2 + self.resultant_vector_y() ** 2) ** 0.5

    # Unit: degree
    def max_spread_direction(self) -> float:
        # Direction of max spread α (alpha) relative to up slope (radian)
        return math.asin(abs(self.resultant_vector_y()) / self.magnitude_of_heading_fire_vector())

    # Unit: ft/min
    def rate_of_spread_heading_fire(self) -> float:
        """Rate of spread of headinging fire in the direction of max spread"""
        return self.rate_of_spread_without_wind_slope() + self.magnitude_of_heading_fire_vector()

    def effective_wind_factor(self) -> float:
        """Effective wind factor in direction of maximum spread"""
        return (self.rate_of_spread_heading_fire() / self.rate_of_spread_without_wind_slope()) - 1

    # Unit: ft/min
    def effective_wind_speed(self) -> float:
        """Effective wind speed in direction of maximum spread"""
        return ((self.effective_wind_factor() * (self.relative_packing_ratio() ** self.eee())) /
                self.ccc()) ** (1/self.bbb())

    def length_to_width_ratio(self):
        """The ellipse length-to-width ratio (z) in the direction of max spread.
        Effective wind speed in mi/h in the formulation of Z."""
        # TODO: check this later
        return 1 + 0.25 * (self.effective_wind_speed * 0.01136)

    def eccentricity(self):
        """The eccentricity of the ellipse (e)"""
        return (((self.length_to_width_ratio() ** 2) - 1) ** 0.5) / self.length_to_width_ratio()

    # Unit: ft/min
    def rate_of_spread_backing_fire(self):
        return self.rate_of_spread_heading_fire() * (1 - self.eccentricity()) / (1 + self.eccentricity())

    # Unit: ft
    def fire_length(self):
        return self.rate_of_spread_heading_fire() + self.rate_of_spread_backing_fire()

    # Unit: ft
    def max_fire_width(self):
        return self.length_to_width_ratio() / self.fire_length()

    # Unit: ft
    def flanking_spread_distance(self):
        return self.max_fire_width() / 2

    # Unit: ft
    def ellipse_dimensions(self):
        fff = self.fire_length() / 2
        ggg = self.rate_of_spread_heading_fire() - fff
        hhh = self.flanking_spread_distance()

        return fff, ggg, hhh

    # Unit: min
    def residence_time(self) -> float:
        return (384 / self.sav_ratio_characteristic()) * self.element_size ** 2

    # Unit: Btu/ft^2
    def heat_per_unit_area(self) -> float:
        return self.reaction_intensity() * self.residence_time()

    # Unit: Btu/ft/s
    def fireline_intensity_heading_fire(self) -> float:
        """Bayram's fire-line intensity. Fire-line intensity for the direction of heading fire"""
        # TODO: Implement this for the directions other than heading fire (check page 91),
        #  it is important because of crown fire implementation.
        return self.heat_per_unit_area() * self.rate_of_spread_heading_fire() / 60

    # Bayram's fire-line intensity for the directions other than heading fire (Btu/ft/s)
    #   TODO: Fill here!!

    # Unit: ft
    def flame_length(self) -> float:
        return 0.45 * self.fireline_intensity_heading_fire() ** 0.46

    # Unit: ft
    def flame_depth(self) -> float:
        return self.rate_of_spread_heading_fire() * self.residence_time()

    # Unit: fraction
    def zex(self, elevation_horizontal_neighbours) -> float:
        """horizontal value for slope steepness and aspect calculation"""
        return -(elevation_horizontal_neighbours[1] - elevation_horizontal_neighbours[0]) / 2 * (self.element_size / 3.281)

    # Unit: fraction
    def zey(self, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        """vertical value for slope steepness and aspect calculation"""
        return -(elevation_vertical_neighbours[1] - elevation_horizontal_neighbours[0]) / 2 * (self.element_size / 3.281)

    # Unit: radian (fraction)
    def calculate_slope_steepness(self, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        """Calculation of slope steepness from elevation"""
        zex2 = self.zex(elevation_horizontal_neighbours) ** 2
        zey2 = self.zey(elevation_horizontal_neighbours, elevation_vertical_neighbours) ** 2
        if zex2 == 0 and zey2 == 0:  # TODO: check this clause later
            return 0
        else:
            return math.sqrt(zex2 + zey2)

    # Unit: degree (azimuth)
    def calculate_aspect(self, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        """Calculation of aspect from elevation"""
        if self.zex(elevation_horizontal_neighbours) == 0 or self.zey(elevation_horizontal_neighbours, elevation_vertical_neighbours) == 0:  # TODO: check this clause later
            aspect = 0
        else:
            aspect = 180 - np.arctan(self.zey(elevation_horizontal_neighbours, elevation_vertical_neighbours) / self.zex(elevation_horizontal_neighbours)) + 90 * (self.zex(elevation_horizontal_neighbours) / abs(self.zex(elevation_horizontal_neighbours)))

        aspect = 180 - aspect  # to fix coordination difference between numpy array and cardinal direction
        aspect += 360 if aspect < 0 else 0  # to get positive degree
        return aspect

    # Unit: degree
    def calculate_upslope_direction(self, elevation_horizontal_neighbours, elevation_vertical_neighbours) -> float:
        """Up-slope direction is the opposite of aspect"""
        upslope_direction = self.calculate_aspect(elevation_horizontal_neighbours, elevation_vertical_neighbours) - 180
        upslope_direction += 360 if upslope_direction < 0 else 0
        return upslope_direction


    def get_element(self, longitude: int, latitude: int) -> Element:

        fuel_model_id: str = self.landscape.fuel_model_distribution[longitude][latitude]
        fuel: Fuel = self.fuel_models[fuel_model_id]
        environment_values: Environment = self.environment[longitude][latitude]
        elevation_horizontal_neighbours = (self.elevation.elevation_distribution[longitude-1][latitude],
                                           self.elevation.elevation_distribution[longitude+1][latitude]
                                           )
        elevation_vertical_neighbours = (self.elevation.elevation_distribution[longitude][latitude-1],
                                         self.elevation.elevation_distribution[longitude][latitude+1]
                                         )


        if fuel is None or fuel_model_id == '0':
            return IncombustibleElement()

        crown_initiation = 0

        return CombustibleElement(
            location=Location(latitude, longitude),
            spread_time=sys.maxsize,
            state=State.FLAMMABLE,
            aspect=environment_values.aspect,
            r_0=self.rate_of_spread_without_wind_slope(),
            r_wind_up_slope=self.rate_of_spread_wind_upslope(),
            residence_time=self.residence_time(),
            heat_per_unit_area=self.heat_per_unit_area(),
            fireline_intensity_heading_fire=self.fireline_intensity_heading_fire(),
            flame_length=self.flame_length(),
            rate_of_spread_headinging_fire=self.rate_of_spread_heading_fire(),
            max_spread_direction=self.max_spread_direction(),
            effective_wind_speed=self.effective_wind_speed(),
            rate_of_spread_backing_fire=self.rate_of_spread_backing_fire(),
            e=self.eccentricity(),
            f=self.ellipse_dimensions()[0],
            g=self.ellipse_dimensions()[1],
            h=self.ellipse_dimensions()[2],
            slope_steepness=self.calculate_slope_steepness(elevation_horizontal_neighbours,
                                                           elevation_vertical_neighbours
                                                           ),
            upslope_direction=self.calculate_upslope_direction(elevation_horizontal_neighbours,
                                                               elevation_vertical_neighbours
                                                               )
        )
