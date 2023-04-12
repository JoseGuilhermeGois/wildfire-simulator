from dataclasses import dataclass


@dataclass
class WindIndex:
    U: list
    V: list
    W: list
    Z_Wnd: list
    ZZ_Wnd: list
    X_Wnd: list
    XX_Wnd: list
    Y_Wnd: list
    YY_Wnd: list


@dataclass
class Environment:
    wind_speed: float
    wind_direction: float
    slope_steepness: float
    aspect: float
