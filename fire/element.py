from dataclasses import dataclass
from enum import Enum
from typing import Protocol

from config.landscape.landscape import Location


class State(Enum):
    NON_FLAMMABLE = 0
    FLAMMABLE = 1
    BURNING = 2
    BURNED = 3


class Element(Protocol):
    state: State


class IncombustibleElement(Element):
    state: State = State.NON_FLAMMABLE


@dataclass
class CombustibleElement(Element):
    state: State
    location: Location
    spread_time: int
    aspect: float
    r_0: float
    r_wind_up_slope: float
    residence_time: float
    heat_per_unit_area: float
    fireline_intensity_head_fire: float
    flame_length: float
    rate_of_spread_heading_fire: float
    max_spread_direction: float
    effective_wind_speed: float
    rate_of_spread_backing_fire: float
    e: float
    f: float
    g: float
    h: float
    upslope_direction: float
    flame_depth: float
    
