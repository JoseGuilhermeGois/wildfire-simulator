from fire import Element, CombustibleElement


def create_file(list_of_ignitions: list[CombustibleElement], ignition_counter, elapsed_time):
    with open("file.out", "w") as outfile:
        outfile.write("V 0.1\n")
        outfile.write(str(elapsed_time) + "    Present time from first ignition\n")
        outfile.write(str(ignition_counter) + "    Total burned cells\n")
        outfile.write("i,   j,    CrdX[m],       CrdY[m], RateofSpread[m/s], ResidenceTime[s], FireDepth[m], "
                      "FireLineInt[kW/m], FlameLength[m], Time[s], Phi[]\n")
        for element in list_of_ignitions:
            outfile.write("{}  {}  {}  {}  {}  {}  {}  {}\n".format(
                element.location.latitude,
                element.location.longitude,
                "{:E}".format(element.rate_of_spread_heading_fire),
                "{:E}".format(element.residence_time),
                "{:E}".format(element.flame_depth),
                "{:E}".format(element.fireline_intensity_head_fire),
                "{:E}".format(element.flame_length * 0.3048),
                "{:E}".format(element.time_of_ignition * 0.3048),
            ))



'''def get_latitude(self, elements: list[list[Element]]):
    return [get_longitude_line_elements(longitude, elements) for longitude in range(self.shape.length)]


def get_longitude_line_elements(self, longitude, elements: list[list[Element]]):
    return [append_to_file(longitude, latitude, elements) for latitude in range(self.shape.width)]


def append_to_file(longitude, latitude, elements):
    element = elements[longitude, latitude]
    if isinstance(element, CombustibleElement) and (element.state.BURNING or element.state.BURNED):
        if element.state.BURNING
            with open("file.out", "a") as f:
                f.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(
                    latitude,
                    longitude,
                    "{:E}".format(element.rate_of_spread_heading_fire),
                    "{:E}".format(element.residence_time),
                    "{:E}".format(element.flame_depth),
                    "{:E}".format(element.fireline_intensity_head_fire),
                    "{:E}".format(element.flame_length * 0.3048),
                    "{:E}".format(element. * 0.3048),
                    "{:E}".format(cell_info.calc_ros()[5] * 3.46146933),
                    int(time_elapsed * 60),
                    "{:E}".format(cell_info.calc_ros()[16]),
                ))
'''