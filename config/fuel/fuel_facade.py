from dataclasses import dataclass
from typing import Protocol

from config.fuel.fuel import Fuel

# Calculations default values
LOAD_CONSTANT = 4.882
MOISTURE_EXTINCTION_DEAD_CONSTANT = 100.0
BED_DEPTH_CONSTANT = 3.28
HEAT_CONTENT_CONSTANT = 2.326
SAV_RATIO_CONSTANT = 30.48


@dataclass
class DefaultsFuel:
    TOTAL_MINERAL_CONTENT_DEFAULT = 0.0555
    EFFECTIVE_MINERAL_CONTENT_DEFAULT = 0.01
    OVEN_DRY_FUEL_LOAD_LIVE_WOOD_DEFAULT = 0.0
    SAV_RATIO_10_H_DEFAULT = 3.576
    SAV_RATIO_100_H_DEFAULT = 0.984
    SAV_RATIO_LIVE_WOOD_DEFAULT = 0.0
    PARTICLE_DENSITY_DEFAULT = 31.96  # 512 divided by 16.018
    MOISTURE_CONTENT_DEFAULT = ((0.1, 0.1, 0.1), (0, 0.80))


# Interface
class FuelFacade(Protocol):

    def get_fuel_models(self, fuel_model_name: str, fuel_model_characteristics: list[float]) -> Fuel:
        ...


class BaseFuelFacade(FuelFacade):

    def __init__(self, defaults_fuel: DefaultsFuel):
        self.defaults = defaults_fuel

    def get_fuel_models(self, fuel_model_name: str, fuel_model_characteristics: list[float]) -> Fuel:
        (
            oven_dry_fuel_load_dead_1h,
            oven_dry_fuel_load_dead_10h,
            oven_dry_fuel_load_dead_100h,
            oven_dry_fuel_load_live_herb,
            fuel_depth_dead,
            fuel_depth_alive,
            surface_area_to_volume_ratio_dead_1h,
            surface_area_to_volume_ratio_alive,
            heat_content_dead,
            heat_content_alive,
            fuel_moisture_dead_1h,
            fuel_moisture_dead_10h,
            fuel_moisture_dead_100h,
            fuel_moisture_alive,
            fuel_moisture_dead_extinction,
            fuel_colour,
            flame_length,
            decay_time,
            initial_rhr_factor
        ) = fuel_model_characteristics

        load = self.to_load(
            oven_dry_fuel_load_dead_1h=oven_dry_fuel_load_dead_1h,
            oven_dry_fuel_load_dead_10h=oven_dry_fuel_load_dead_10h,
            oven_dry_fuel_load_dead_100h=oven_dry_fuel_load_dead_100h,
            oven_dry_fuel_load_live_herb=oven_dry_fuel_load_live_herb)

        sav_ratio = self.to_sav_ratio(
            sav_ratio_1_h=surface_area_to_volume_ratio_dead_1h,
            sav_ratio_live_herb=surface_area_to_volume_ratio_alive)

        heat_content = self.to_heat_content(
            heat_content_dead=heat_content_dead,
            heat_content_alive=heat_content_alive)

        bed_depth = self.to_bed_depth(fuel_depth_dead)

        moisture_extinction_dead = self.to_moisture_extinction_dead(fuel_moisture_dead_extinction)

        return Fuel(
            name=fuel_model_name,
            load=load,
            sav_ratio=sav_ratio,
            heat_content=heat_content,
            bed_depth=bed_depth,
            moisture_extinction_dead=moisture_extinction_dead,
            particle_density=self.defaults.PARTICLE_DENSITY_DEFAULT,
            total_mineral_content=self.defaults.TOTAL_MINERAL_CONTENT_DEFAULT,
            effective_mineral_content=self.defaults.EFFECTIVE_MINERAL_CONTENT_DEFAULT,
            moisture_content=self.defaults.MOISTURE_CONTENT_DEFAULT
        )

    def to_load(self,
                oven_dry_fuel_load_dead_1h: float,
                oven_dry_fuel_load_dead_10h: float,
                oven_dry_fuel_load_dead_100h: float,
                oven_dry_fuel_load_live_herb: float) -> tuple:

        raw_fuel_load_dead = (oven_dry_fuel_load_dead_1h, oven_dry_fuel_load_dead_10h, oven_dry_fuel_load_dead_100h)
        fuel_load_dead = tuple(map(lambda x: x / LOAD_CONSTANT, raw_fuel_load_dead))

        raw_fuel_load_live = (oven_dry_fuel_load_live_herb, self.defaults.OVEN_DRY_FUEL_LOAD_LIVE_WOOD_DEFAULT)
        fuel_load_live = tuple(map(lambda x: x / LOAD_CONSTANT, raw_fuel_load_live))

        return tuple([fuel_load_dead, fuel_load_live])

    def to_sav_ratio(self, sav_ratio_1_h: float, sav_ratio_live_herb: float):
        raw_sav_ratio_dead = (sav_ratio_1_h, self.defaults.SAV_RATIO_10_H_DEFAULT, self.defaults.SAV_RATIO_100_H_DEFAULT)
        sav_ratio_dead = tuple(map(lambda x: x * SAV_RATIO_CONSTANT, raw_sav_ratio_dead))

        raw_sav_ratio_live = (sav_ratio_live_herb, self.defaults.SAV_RATIO_LIVE_WOOD_DEFAULT)
        sav_ratio_live = tuple(map(lambda x: x * SAV_RATIO_CONSTANT, raw_sav_ratio_live))

        return tuple([sav_ratio_dead, sav_ratio_live])

    @staticmethod
    def to_heat_content(heat_content_dead: float, heat_content_alive: float):
        raw_heat_content = (heat_content_dead, heat_content_alive)
        heat_content = tuple(map(lambda x: x / HEAT_CONTENT_CONSTANT, raw_heat_content))

        return heat_content

    @staticmethod
    def to_bed_depth(bed_depth: float) -> float:
        return float(bed_depth * BED_DEPTH_CONSTANT)

    @staticmethod
    def to_moisture_extinction_dead(moisture_extinction_dead: float) -> float:
        return float(moisture_extinction_dead / MOISTURE_EXTINCTION_DEAD_CONSTANT)
