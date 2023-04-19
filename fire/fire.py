from config import Landscape
from fire import Element, State, CombustibleElement, Frame, visualize
from .spread import Neighbour, SpreadTimeCalculatorStrategyFactory


class Fire:
    def __init__(self, landscape: Landscape, terrain_topography: list[list[Element]], iterations: int):
        self.shape = landscape.shape
        self.element_size: int = landscape.element_size
        self.terrain_topography: list[list[Element]] = terrain_topography
        self.iterations = iterations

    def start(self):
        frames = [self.frame() for _ in range(self.iterations)]
        visualize(frames)

    def frame(self):
        lowest_spread_time = self.update_spread_and_get_lowest()
        states = self.step_forward_and_return_new_states(time_interval=lowest_spread_time)
        return Frame(data=states, interval=lowest_spread_time)

    """ Calculate spread! """

    def update_spread_and_get_lowest(self) -> float:
        return min(
            [self.update_longitude_line_spread_times_and_get_lowest(longitude)
             for longitude in range(self.shape.length)]
        )

    def update_longitude_line_spread_times_and_get_lowest(self, longitude: int) -> float:
        times = [self.update_coordinate_spread_time(latitude, longitude) for latitude in range(self.shape.width)]
        return min(filter(lambda item: item is not None, times))

    def update_coordinate_spread_time(self, latitude: int, longitude: int) -> float | None:
        element = self.terrain_topography[latitude][longitude]
        if isinstance(element, CombustibleElement) and element.state == State.FLAMMABLE:
            spread_time_calculator_strategy_factory = SpreadTimeCalculatorStrategyFactory(element, self.element_size)
            lowest_spread_time = element.spread_time
            for neighbour in Neighbour:
                spread_time_calculator_strategy = spread_time_calculator_strategy_factory.create(neighbour)
                for neighbour_latitude_offset, neighbour_longitude_offset in neighbour.value:
                    neighbour_latitude = longitude + neighbour_latitude_offset
                    neighbour_longitude = latitude + neighbour_longitude_offset
                    if not self.is_valid_coordinate(neighbour_latitude, neighbour_longitude):
                        continue

                    neighbour_element = self.terrain_topography[neighbour_latitude][neighbour_longitude]
                    if neighbour_element.state == State.BURNING:
                        spread_time = spread_time_calculator_strategy.calculate(neighbour_element)
                        if spread_time is not None and spread_time < lowest_spread_time:
                            lowest_spread_time = spread_time
            element.spread_time = lowest_spread_time
            return element.spread_time

    """ Steps to go forward! To the next step. """

    def step_forward_and_return_new_states(self, time_interval: float) -> list[list[float]]:
        return [self.update_latitude_line_elements(latitude, time_interval) for latitude in range(self.shape.width)]

    def update_latitude_line_elements(self, latitude: int, time: float) -> list[float]:
        return [self.update_coordinate_element(latitude, longitude, time) for longitude in range(self.shape.length)]

    def update_coordinate_element(self, latitude: int, longitude: int, time_interval: float):
        element = self.terrain_topography[latitude][longitude]
        if isinstance(element, CombustibleElement):
            return self.update_combustible_element(element, time_interval)
        return element.state.value

    def update_combustible_element(self, element: CombustibleElement, time_interval: float):
        if element.state == State.FLAMMABLE:
            element.spread_time -= time_interval
            if element.spread_time <= 0:
                element.state = State.BURNING

        elif element.state == State.BURNING:
            element.residence_time -= time_interval
            if element.residence_time <= 0:
                element.state = State.BURNED

        return element.state.value

    def is_valid_coordinate(self, latitude: int, longitude: int) -> bool:
        return latitude in range(self.shape.width) and longitude in range(self.shape.length)
