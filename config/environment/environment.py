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


class Environment:
