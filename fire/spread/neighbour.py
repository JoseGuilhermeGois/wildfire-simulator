from enum import Enum


class Neighbour(Enum):
    ADJACENT = ((-1, 0), (0, -1), (0, 1), (1, 0))
    DIAGONAL = ((-1, -1), (-1, 1), (1, -1), (1, 1))
    DISTANT = ((-1, 2), (1, 2), (2, 1), (2, -1), (1, -2), (-1, -2), (-2, -1), (-2, 1))
