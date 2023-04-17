import math
from abc import ABC, abstractmethod
from fire.element import Element
from fire.spread import Neighbour, get_rate_of_spread


class SpreadTimeCalculatorStrategyFactory:

    def __init__(self, central_element: Element, element_size: float):
        self.central_element = central_element
        self.element_size = element_size

    def create(self, neighbour_type: Neighbour):
        if neighbour_type == Neighbour.ADJACENT:
            return AdjacentNeighbourSpreadTimeCalculator(self.central_element, self.element_size)
        elif neighbour_type == Neighbour.DIAGONAL:
            return DiagonalNeighbourSpreadTimeCalculator(self.central_element, self.element_size)
        elif neighbour_type == Neighbour.DISTANT:
            return DistantNeighbourSpreadTimeCalculator(self.central_element, self.element_size)
        else:
            raise NotImplementedError


class SpreadTimeCalculatorStrategy(ABC):

    def __init__(self, central_element: Element, element_size: float):
        self.central_element = central_element
        self.element_size = element_size

    @abstractmethod
    def distance(self) -> float:
        raise NotImplementedError

    def calculate(self, neighbour_element: Element) -> float:
        rate_of_spread = get_rate_of_spread(element1=self.central_element, element2=neighbour_element)
        if not rate_of_spread:
            return self.distance() / rate_of_spread

        # return self.central_element.reset_ignition_time()


class AdjacentNeighbourSpreadTimeCalculator(SpreadTimeCalculatorStrategy):

    def __init__(self, central_element: Element, element_size: float):
        super().__init__(central_element, element_size)

    def distance(self) -> float:
        return self.element_size


class DiagonalNeighbourSpreadTimeCalculator(SpreadTimeCalculatorStrategy):

    def __init__(self, central_element: Element, element_size: float):
        super().__init__(central_element, element_size)

    def distance(self) -> float:
        return math.sqrt(2) * self.element_size


class DistantNeighbourSpreadTimeCalculator(SpreadTimeCalculatorStrategy):

    def __init__(self, central_element: Element, element_size: float):
        super().__init__(central_element, element_size)

    def distance(self) -> float:
        return math.sqrt((self.element_size ** 2) + ((2 * self.element_size) ** 2))
