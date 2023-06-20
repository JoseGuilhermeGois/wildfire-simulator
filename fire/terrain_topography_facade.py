import math
import sys
from typing import Protocol

import numpy as np

from config.reader.environment import EnvironmentValues
from config.reader.environment.environment import WindIndex
from config.reader.fuel import Fuel
from fire.element import CombustibleElement, IncombustibleElement, Element, State
from config.reader.landscape import Landscape, Location
from config.reader.elevation import Elevation


class TerrainTopographyFacade(Protocol):

    def get_element(self, longitude: int, latitude: int) -> Element:
        ...


class BaseTerrainTopographyFacade(TerrainTopographyFacade):

    def __init__(self, landscape: Landscape, fuel_models: dict[str, Fuel], environment: WindIndex,
                 elevation: Elevation):
        self.landscape: Landscape = landscape
        self.elevation: Elevation = elevation
        self.fuel_models: dict[str, Fuel] = fuel_models
        self.environment: WindIndex = environment

    def get_element(self, y: int, x: int) -> Element:

        fuel_model_id: str = self.landscape.fuel_model_distribution[y][x]
        fuel: Fuel = self.fuel_models[fuel_model_id]

        if fuel is None or fuel_model_id == '0':
            return IncombustibleElement()

        # Weighting factors calculations
        # Mean total surface area per unit fuel cell of each size class within each category
        area_ij = [[0] * 3, [0] * 2]
        for i in range(len(area_ij)):
            for j in range(len(area_ij[i])):
                area_ij[i][j] = fuel.sav_ratio[i][j] * fuel.load[i][j] / fuel.particle_density

        # Mean total surface area of the live and dead categories
        area_i = (sum(area_ij[0]), sum(area_ij[1]))
        area_t = sum(area_i)

        # Weighting factor for characteristic dead and live heat content, effective mineral content,
        # moisture content and surface-area-to-volume ratio
        # factor_ij: [[dead_1h, dead_10h, dead_100h], [live_herb, live_woody]]
        factor_ij = [[0] * 3, [0] * 2]
        for i in range(len(factor_ij)):
            for j in range(len(factor_ij[i])):
                factor_ij[i][j] = area_ij[i][j] / area_i[i] if area_i[i] > 0 else 0

        # Weighting factor for characteristic fuel bed surface-area-to-volume ratio
        # f_i: (dead, live)
        factor_i = (area_i[0] / area_t, area_i[1] / area_t)

        # Characteristic values
        # net_load: [[dead_1h, dead_10h, dead_100h], [live_herb, live_woody]]
        net_load = [[0] * 3, [0] * 2]
        for i in range(len(net_load)):
            for j in range(len(net_load[i])):
                net_load[i][j] = fuel.load[i][j] * (1 - fuel.total_mineral_content)

        # omg0d in Firestation
        net_load_dead = \
            factor_ij[0][0] * net_load[0][0] + factor_ij[0][1] * net_load[0][1] + factor_ij[0][2] * net_load[0][2]

        net_load_live = net_load[1][0] + net_load[1][1]

        heat_content_i = (sum(factor_ij[0]) * fuel.heat_content[0], sum(factor_ij[1]) * fuel.heat_content[1])

        if sum(factor_ij[1]) > 0:
            effective_mineral_content_i = (sum(factor_ij[0]) * fuel.effective_mineral_content,
                                           sum(factor_ij[1]) * fuel.effective_mineral_content)

            mineral_damping_coff_i = (min(0.174 * effective_mineral_content_i[0] ** (-0.19), 1),
                                      min(0.174 * effective_mineral_content_i[1] ** (-0.19), 1))
        else:
            effective_mineral_content_i = (sum(factor_ij[0]) * fuel.effective_mineral_content, 0)
            mineral_damping_coff_i = (min(0.174 * effective_mineral_content_i[0] ** (-0.19), 1), 0)

        # aamfd in Firestation
        moisture_content_i = (factor_ij[0][0] * fuel.moisture_content[0][0] +
                              factor_ij[0][1] * fuel.moisture_content[0][1] +
                              factor_ij[0][2] * fuel.moisture_content[0][2],
                              factor_ij[1][0] * fuel.moisture_content[1][0] +
                              factor_ij[1][1] * fuel.moisture_content[1][1])

        # Moisture Extinction
        dead_to_live_load_ratio_num = 0
        fine_dead_fuel_moisture_den = 0
        for j in range(len(fuel.sav_ratio[0])):
            if fuel.sav_ratio[0][j] > 0:
                dead_to_live_load_ratio_num += fuel.load[0][j] * math.exp(-138 / (fuel.sav_ratio[0][j]))
                fine_dead_fuel_moisture_den += fuel.load[0][j] * math.exp(-138 / (fuel.sav_ratio[0][j]))

        dead_to_live_load_ratio_den = 0
        for j in range(len(fuel.sav_ratio[1])):
            if fuel.sav_ratio[1][j] > 0:
                dead_to_live_load_ratio_den += fuel.load[1][j] * math.exp(-500 / (fuel.sav_ratio[1][j]))

        if dead_to_live_load_ratio_den > 0:
            dead_to_live_load_ratio = dead_to_live_load_ratio_num / dead_to_live_load_ratio_den
        else:
            dead_to_live_load_ratio = 0

        fine_dead_fuel_moisture_num = 0
        for j in range(len(fuel.sav_ratio[0])):
            if fuel.sav_ratio[0][j] > 0:
                fine_dead_fuel_moisture_num += fuel.moisture_content[0][j] * fuel.load[0][j] * \
                                            math.exp(-138 / (fuel.sav_ratio[0][j]))

        fine_dead_fuel_moisture = fine_dead_fuel_moisture_num / fine_dead_fuel_moisture_den if \
            fine_dead_fuel_moisture_den > 0 else 0

        moisture_extinction_live = max(2.9 * dead_to_live_load_ratio *
                                       (1 - (fine_dead_fuel_moisture / fuel.moisture_extinction_dead)) - 0.226,
                                       fuel.moisture_extinction_dead)

        tau_m = [min(moisture_content_i[0] / fuel.moisture_extinction_dead, 1),
                 min(moisture_content_i[1] / moisture_extinction_live, 1)]

        moisture_damping_coff_i = [0] * 2
        for i in range(0, 2):
            moisture_damping_coff_i[i] = 1 - 2.59 * tau_m[i] + 5.11 * (tau_m[i] ** 2) - 3.52 * (tau_m[i] ** 3)

        # Sigmad(n) in Firestation
        sav_ratio_i = (factor_ij[0][0] * fuel.sav_ratio[0][0] + factor_ij[0][1] * fuel.sav_ratio[0][1] +
                       factor_ij[0][2] * fuel.sav_ratio[0][2], factor_ij[1][0] * fuel.sav_ratio[1][0] +
                       factor_ij[1][1] * fuel.sav_ratio[1][1])

        sav_ratio_characteristic = factor_i[0] * sav_ratio_i[0] + factor_i[1] * sav_ratio_i[1]

        mean_bulk_density = (sum(fuel.load[0]) + sum(fuel.load[1])) / fuel.bed_depth

        mean_packing_ratio = mean_bulk_density / fuel.particle_density

        optimum_packing_ratio = 3.348 * sav_ratio_characteristic ** (-0.8189)

        relative_packing_ratio = mean_packing_ratio / optimum_packing_ratio

        # Heat Source
        propagating_flux_ratio = (192 + 0.2595 * sav_ratio_characteristic) ** (-1) * math.exp(
            (0.792 + 0.681 * sav_ratio_characteristic ** 0.5) * (mean_packing_ratio + 0.1))

        aaa = 133 * sav_ratio_characteristic ** (-0.7913)

        maximum_reaction_velocity = \
            sav_ratio_characteristic ** 1.5 * (495 + 0.0594 * sav_ratio_characteristic ** 1.5) ** (-1)

        optimum_reaction_velocity = \
            maximum_reaction_velocity * (relative_packing_ratio ** aaa) * math.exp(aaa * (1 - relative_packing_ratio))

        reaction_intensity_dead = net_load_dead * heat_content_i[0] * moisture_damping_coff_i[0] * \
                                  mineral_damping_coff_i[0]
        reaction_intensity_live = net_load_live * heat_content_i[1] * moisture_damping_coff_i[1] * \
                                  mineral_damping_coff_i[1]

        reaction_intensity = optimum_reaction_velocity * (reaction_intensity_dead + reaction_intensity_live)

        heat_of_pre_ignition_ij = [[0] * 3, [0] * 2]
        for i in range(2):
            for j in range(len(heat_of_pre_ignition_ij[i])):
                heat_of_pre_ignition_ij[i][j] = 250 + 1116 * fuel.moisture_content[i][j]

        heat_number_tmp = [0] * 2
        for i in range(len(heat_number_tmp)):
            heat_number_tmp[i] = \
                factor_i[i] * sum([a * (math.exp(-138 / b) if b > 0 else 0) * c for a, b, c in zip(
                    factor_ij[i], fuel.sav_ratio[i], heat_of_pre_ignition_ij[i])])

        heat_sink = mean_bulk_density * sum(heat_number_tmp)

        ccc = 7.47 * math.exp(-0.133 * (sav_ratio_characteristic ** 0.55))
        bbb = 0.02526 * (sav_ratio_characteristic ** 0.54)
        eee = 0.715 * math.exp(-0.000359 * sav_ratio_characteristic)

        '''///////////////////////////// CALCULATE WIND //////////////////////////////////////////'''
        crd_x_node = self.landscape.location.real_latitude + (
                x - 1) * self.landscape.element_size + 0.5 * self.landscape.element_size
        crd_y_node = self.landscape.location.real_longitude + (
                y - 1) * self.landscape.element_size + 0.5 * self.landscape.element_size

        vel_u, vel_v, vel_w = self.convert_crd_to_wind_mesh(crd_x_node, crd_y_node, fuel)
        wind_angle = self.calc_angle(vel_u, vel_v)
        vel_wind = math.sqrt(vel_u * vel_u + vel_v * vel_v + 0 * vel_w * vel_w)
        valor_limit = 0.0004023 * reaction_intensity
        # It is now recommended that a wind limit not be imposed (Andrews et al. 2013) (including Rothermel).
        # wind_limit = min(self.wind_speed, 96.8 * reaction_intensity ** 1 / 3)
        if vel_wind > valor_limit:
            vel_wind = valor_limit
        vel_wind_ft_min = vel_wind * 196.82
        '''///////////////////////////////////////////////////////////////////////////////////////'''

        pi_w = ccc * (vel_wind_ft_min ** bbb) * (relative_packing_ratio ** (-eee))

        angle = math.pi * self.elevation.slope_steepness_list[y+1][x+1] / 1800

        pi_s = 5.275 * (mean_packing_ratio ** (-0.3)) * (math.tan(angle) ** 2)

        heat_source = reaction_intensity * propagating_flux_ratio * (1 + pi_w + pi_s)

        # Rate of spread with no wind and slope (ft/min)
        r_0 = reaction_intensity * propagating_flux_ratio / heat_sink

        # Rate of spread if wind is up slope (ft/min)
        r_wind_up_slope = heat_source / heat_sink
        r_wind_up_slope = round(r_wind_up_slope, 6)

        # ----------------------------------------------------------------------------------------
        d_s = r_0 * pi_s  # * t (elapsed time)
        d_w = r_0 * pi_w  # * t (elapsed time)
        ''' 
        The wind direction used in this model is the opposite of the given wind direction. 
        Because the given wind direction is defined as the direction the wind is coming from.
        For this model, we need the direction the wind goes.
        '''
        wind_direction = wind_angle
        # wind_direction += 360 if wind_direction < 0 else 0
        '''
        The upslope direction is the opposite of aspect. 
        '''
        # upslope_direction = self.elevation.up_slope_direction_list[y][x]  # - 180
        # upslope_direction += 360 if upslope_direction < 0 else 0
        upslope_direction = 0.1 * math.pi * self.elevation.up_slope_direction_list[y+1][x+1] / 180

        phi_x = pi_w * math.cos(wind_angle) + pi_s * math.cos(upslope_direction)
        phi_y = pi_w * math.sin(wind_angle) + pi_s * math.sin(upslope_direction)
        phi = math.sqrt(phi_x * phi_x + phi_y * phi_y)
        angle_velocity_equiv = math.atan2(phi_y, phi_x)
        prop_dir = angle_velocity_equiv * 180 / math.pi
        # effective_wind_speed = (((phi * (relative_packing_ratio ** eee)) / ccc) ** (1 / bbb))
        # vel_equiv_mih = effective_wind_speed * 0.01136
        # vel_equiv_kmh = effective_wind_speed * 0.01829

        # Wind direction ω (omega) relative to upslope (degree)
        omega = wind_direction - upslope_direction
        omega += 360 if omega < 0 else 0

        xx = d_s + (d_w * math.cos(math.radians(omega)))
        yy = d_w * math.sin(math.radians(omega))

        d_h = (xx ** 2 + yy ** 2) ** 0.5

        # Direction of max spread α (alpha) relative to up slope (radian)
        max_spread_direction = math.asin(abs(yy) / d_h) if d_h > 0 else 0

        # Rate of spread of heading fire in the direction of max spread (ft/min)
        r_h = r_0 + d_h  # * t (elapsed time)
        r_h = round(r_h, 6)

        effective_wind_factor = (r_h / r_0) - 1

        # Effective wind speed (ft/min)
        effective_wind_speed = (((effective_wind_factor * (relative_packing_ratio ** eee)) / ccc) ** (1 / bbb)) \
            if ccc > 0 and bbb > 0 else 0

        # The ellipse length-to-width ratio in the direction of max spread.
        # Effective wind speed in mi/h in the formulation of Z. # TODO: check this later
        zzz = 1 + 0.25 * (effective_wind_speed * 0.01136)

        # The eccentricity of the ellipse
        e = (((zzz ** 2) - 1) ** 0.5) / zzz

        # Rate of spread of backing fire (ft/min)
        r_b = r_h * (1 - e) / (1 + e)

        # Fire length
        lll = r_h + r_b
        # Max fire width
        www = zzz / lll if lll != 0 else 0
        # Flanking spread distance
        d_f = www / 2

        # Ellipse dimensions
        f = lll / 2
        g = r_h - f
        h = d_f

        # Residence time (min)
        tr = (384 / sav_ratio_characteristic) * (self.landscape.element_size ** 2)
        tr = round(tr, 3)

        # Heat per unit area (Btu/ft^2)
        ha = reaction_intensity * tr

        # Bayram's fireline intensity. Fireline intensity for the direction of head fire (Btu/ft/s)
        #   TODO: Implement this for the directions other than head fire (check page 91).
        #    Important because of crown fire implementation.
        fire_line_intensity_head_fire = ha * r_h / 60

        # Bayram's fireline intensity for the directions other than head fire (Btu/ft/s)
        #   TODO: Fill here!!

        # Flame length (ft)
        flame_length = 0.45 * fire_line_intensity_head_fire ** 0.46
        flame_depth = r_h * tr

        ''' IMPLEMENTATION OF CROWN FIRE '''
        '''# Crown Fire Initiation
        i_0 = (0.010 * self.cbh * (460 + 25.9 * self.fmc)) ** (3 / 2)
        # Crown Bulk Density
        CFL = 1.19
        if self.stand_height == 0 or self.stand_height <= self.cbh:
            CBD = 0
        else:
            CBD = CFL / (self.stand_height - self.cbh)
        # Wind 10 meter open-wind
        U_10 = 20
        # Estimate Fine Fuel Moisture
        EFFM = 9

        if self.stand_height > 0:
            if (fire_line_intensity_head_fire * 3.461469) >= i_0:
                crown_initiation = 1

                # Cruz, Alexandre and Wakimoto (2005)
                # Threshold for active spread fire
                if CBD != 0:
                    CAC = 3 / CBD
                else:
                    CAC = 0

                # Rate of crown fire spread
                if CAC >= 1.0:  # Active Crown Fire
                    CROS_a = 11.02 * U_10**0.9 * CBD**0.19 * np.exp(-0.17 * EFFM)
                    if (3.28084 * CROS_a) > r_h:
                        r_h = CROS_a * 3.28084
                else:  # Passive Crown Fire
                    CROS_p = (11.02 * U_10**0.9 * CBD**0.19 * np.exp(-0.17 * EFFM)) * np.exp(-CAC)
                    if (CROS_p * 3.28084) > r_h:
                        r_h = CROS_p * 3.28084
            else:
                crown_initiation = 0
        else:
            crown_initiation = 0'''
        crown_initiation = 0

        return CombustibleElement(
            location=Location(crd_x_node, crd_y_node),
            latitude=x,
            longitude=y,
            spread_time=sys.maxsize,
            state=State.FLAMMABLE,
            slope_steepness=self.elevation.slope_steepness_list[y+1][x+1],
            r_0=r_0,
            r_wind_up_slope=r_wind_up_slope,
            fixed_residence_time=tr,
            residence_time=tr,
            heat_per_unit_area=ha,
            fireline_intensity_head_fire=fire_line_intensity_head_fire,
            flame_length=flame_length,
            rate_of_spread_heading_fire=r_h,
            max_spread_direction=max_spread_direction,
            effective_wind_speed=effective_wind_speed,
            rate_of_spread_backing_fire=r_b,
            e=e,
            f=f,
            g=g,
            h=h,
            upslope_direction=upslope_direction,
            flame_depth=flame_depth,
            time_of_ignition=None,
            ros=None
        )

    def convert_crd_to_wind_mesh(self, crd_x_node, crd_y_node, fuel):
        iv = int((crd_x_node - self.environment.X_IniWnd) / self.environment.DxWnd + 0.5)
        jv = int((crd_y_node - self.environment.Y_IniWnd) / self.environment.DyWnd + 0.5)
        if iv < 0 or iv >= self.environment.NiWnd or jv < 0 or jv >= self.environment.NjWnd:
            vel_u = 0
            vel_v = 0
            vel_w = 0
            return vel_u, vel_v, vel_w

        xv1 = self.environment.XX_Wnd[iv]
        xv2 = self.environment.XX_Wnd[iv + 1]
        yv1 = self.environment.YY_Wnd[jv]
        yv2 = self.environment.YY_Wnd[jv + 1]

        height = float(self.environment.ZZ_Wnd[iv][jv][1])
        flngt = fuel.flame_length
        veglngt = fuel.bed_depth

        vel_u = self.get_velocity_component(self.environment.U, iv, jv, xv1, xv2, yv1, yv2, crd_x_node, crd_y_node, height,
                                            flngt, veglngt)
        vel_v = self.get_velocity_component(self.environment.V, iv, jv, xv1, xv2, yv1, yv2, crd_x_node, crd_y_node, height,
                                            flngt, veglngt)
        vel_w = self.get_velocity_component(self.environment.W, iv, jv, xv1, xv2, yv1, yv2, crd_x_node, crd_y_node, height,
                                            flngt, veglngt)

        return vel_u, vel_v, vel_w

    def get_velocity_component(self, wind_component, iv, jv, xv1, xv2, yv1, yv2, crd_x_node, crd_y_node, height, flngt,
                               veglngt):

        velfga = self.interpol(yv1, wind_component[iv][jv][1], yv2, wind_component[iv][jv + 1][1], crd_y_node)
        velfgb = self.interpol(yv1, wind_component[iv + 1][jv][1], yv2, wind_component[iv + 1][jv + 1][1], crd_y_node)
        result = self.interpol(xv1, velfga, xv2, velfgb, crd_x_node)
        if ((height + 0.36 * veglngt) / (0.13 * veglngt)) <= 0:
            vel_comp = 0
        else:
            vel_comp = (result * (1 + 0.36 * veglngt / flngt) / (math.log((height + 0.36 * veglngt) / (0.13 * veglngt))))*(
                math.log((flngt / veglngt + 0.36) / 0.13) - 1)

        return vel_comp

    def interpol(self, x1, y1, x2, y2, x):
        return y1 + (y2 - y1) / (0.0000000001 + x2 - x1) * (x - x1)

    @staticmethod
    def calc_angle(number1, number_2):
        if number1 > 0 and number_2 >= 0:
            return math.atan(number_2 / number1)
        elif number1 < 0 <= number_2:
            return math.atan(number_2 / number1) + 3.14159
        elif number1 < 0 and number_2 <= 0:
            return math.atan(number_2 / number1) + 3.14159
        elif number1 > 0 > number_2:
            return math.atan(number_2 / number1) + 2 * 3.14159
        elif number1 == 0 and number_2 > 0:
            return 0.5 * 3.14159
        elif number1 == 0 and number_2 < 0:
            return -0.5 * 3.14159
        elif number1 == 0 and number_2 == 0:
            return 0
