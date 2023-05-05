from dataclasses import dataclass
from typing import Protocol


class Environment(Protocol):
    pass


@dataclass
class WindIndex(Environment):
    type: str
    U: list
    V: list
    W: list
    Z_Wnd: list
    ZZ_Wnd: list
    X_Wnd: list
    XX_Wnd: list
    Y_Wnd: list
    YY_Wnd: list
    DxWnd: float
    DyWnd: float
    X_IniWnd: float
    Y_IniWnd: float
    X_FinWnd: float
    Y_FinWnd: float
    NiWnd: int
    NjWnd: int
    NkWnd: int
    AuxSingle: float


@dataclass
class WindDefault(Environment):
    wind_speed: float
    wind_direction: float
    slope_steepness: float
    aspect: float
