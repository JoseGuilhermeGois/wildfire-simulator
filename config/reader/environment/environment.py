from dataclasses import dataclass


@dataclass
class WindIndex:
    type: str
    U: list[list[list[float]]]
    V: list[list[list[float]]]
    W: list[list[list[float]]]
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
class WindDefault:
    wind_speed: float
    wind_direction: float
    slope_steepness: float
    aspect: float


@dataclass
class EnvironmentValues:
    slope_steepness: float
    aspect: float
