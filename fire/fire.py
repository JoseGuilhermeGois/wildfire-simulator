from simulator.field import FieldState
from simulator.topography.grid import get_state_grid, get_fuel_grid, create_grid
from spread import Neighbour

TOTAL_TIME = 100

# ---------------- IMPORTS (OWN FILES) -----------------

from simulator.topography.topography import *

grid = create_grid([(57, 10)])

# ---------------- FIRE SPREAD MODEL -------------------

results = []  # For simulation

print("States at the beginning:")
print(get_state_grid(grid))
print()
print("Fuel Map:")
print(get_fuel_grid(grid))

elapsed_time_collection = []
elapsed_time = 0.0
elapsed_time_collection.append(elapsed_time)

for time in range(0, TOTAL_TIME):

    ignition_times_of_all_cells = []

    calculate_current_ignition_times()

    print()
    print("Ignition times of all cells:")
    print(ignition_times_of_all_cells)

    if min(ignition_times_of_all_cells) == sys.maxsize:
        shortest_ignition_time = 0
    else:
        shortest_ignition_time = min(ignition_times_of_all_cells)

    print()
    print("Shortest ignition time of this time step(min):")
    print(shortest_ignition_time)

    update_grid()

    elapsed_time += shortest_ignition_time
    print()
    print("Elapsed time(min):")
    print(elapsed_time)
    elapsed_time_collection.append(round(elapsed_time))

    print()
    print("States at the end of time: " + str(time))
    print(get_state_grid(grid))
    print()

    results.append(get_state_grid(grid))  # For simulation

print("Fuel Map:")
print(get_fuel_grid(grid))


def calculate_current_ignition_times():

    # Ignore edges in calculation
    for y in range(2, GRID_HEIGHT - 2):
        for x in range(2, GRID_WIDTH - 2):

            current_cell = grid[y][x]

            print()
            print("----rate of spread heading fire(ft/min):----")
            print(current_cell.rate_of_spread_heading_fire)

            if current_cell == FieldState.BURNED:
                current_cell.reset_ignition_time()

            elif current_cell.state == FieldState.FLAMMABLE:

                current_cell.has_burning_neighbor = False

                # List of all ignition times of current cell
                ignition_times_of_current_cell = []

                ignitionTimeCalculatorStrategyFactory = IgnitionTimeCalculatorStrategyFactory(current_cell)
                for neighbour in Neighbour:
                    ignition_time_calculator = ignitionTimeCalculatorStrategyFactory.create(neighbour)
                    for neighbour_row, neighbour_column in neighbour.value:
                        neighbour_cell = grid[y + neighbour_row][x + neighbour_column]
                        if neighbour_cell.state == FieldState.BURNING:
                            ignition_time = ignition_time_calculator.calculate(neighbour_cell)
                            ignition_times_of_current_cell.append(ignition_time)

                if ignition_times_of_current_cell:
                    ignition_times_of_current_cell.append(current_cell.ignition_time)
                    print("Ignition times of current cell:")
                    print(ignition_times_of_current_cell)

                    current_cell.ignition_time = min(ignition_times_of_current_cell)
                    ignition_times_of_all_cells.append(current_cell.ignition_time)
                #
                print("Ignition time of current cell(min):")
                print(current_cell.ignition_time)


def update_grid():
    for y in range(2, GRID_HEIGHT - 2):
        for x in range(2, GRID_WIDTH - 2):

            current_cell = grid[y][x]

            # Update ignition time
            if current_cell.ignition_time is not None:
                current_cell.ignition_time = current_cell.ignition_time - shortest_ignition_time

            # Update state
            if current_cell.state == FieldState.BURNING:

                if current_cell.residence_time <= shortest_ignition_time:
                    current_cell.set_state(FieldState.BURNED)
                    current_cell.residence_time = None
                    current_cell.reset_ignition_time()

                elif current_cell.residence_time > shortest_ignition_time:

                    current_cell.residence_time -= shortest_ignition_time

            elif current_cell.state == FieldState.FLAMMABLE and current_cell.ignition_time <= 0 and current_cell.ignition_time != None:
                current_cell.set_state(FieldState.BURNING)

            # print()
            # print('ignition_time_after_reducing_shortest_ignition_time:')
            # print(current_cell.ignition_time)
