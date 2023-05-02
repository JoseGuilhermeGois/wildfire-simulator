from dataclasses import dataclass


@dataclass
class Location:
    real_latitude: float
    real_longitude: float


@dataclass
class Shape:
    width: int
    length: int


@dataclass
class Landscape:
    shape: Shape
    location: Location
    element_size: int
    fuel_model_distribution: list[list[str]]
