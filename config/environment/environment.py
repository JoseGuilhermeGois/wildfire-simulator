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
class NodesNumber:
    NiWnd: int
    NjWnd: int
    NkWnd: int
    AuxSingle: float


@dataclass
class WindMeshSpacing:
    DxWnd: float
    DyWnd: float


@dataclass
class Environment:
    X_IniWnd: float
    Y_IniWnd: float
    X_FinWnd: float
    Y_FinWnd: float
