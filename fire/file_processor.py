from fire import CombustibleElement


def create_file(list_of_ignitions: list[CombustibleElement], ignition_counter, elapsed_time):
    with open("C:\\Users\\guigo\\Software\\IMfireProject\\output\\fire.out", "w") as outfile:
        outfile.write("V 0.1\n")
        outfile.write(str(round(elapsed_time * 60, 1)) + "    Present time from first ignition in seconds\n")
        outfile.write(str(ignition_counter) + "    Total burned cells\n")
        outfile.write("i,   j,    CrdX[m],       CrdY[m], RateofSpread[m/s], ResidenceTime[s], FireDepth[m], "
                      "FireLineInt[kW/m], FlameLength[m], Time[s], Phi[]\n")
        for element in list_of_ignitions:
            outfile.write("{}  {}  {}  {}  {}  {}  {}  {}  {}  {}\n".format(
                element.latitude,
                element.longitude,
                element.location.real_latitude,
                element.location.real_longitude,
                "{:E}".format(element.rate_of_spread_heading_fire * 0.00508),
                "{:E}".format(element.fixed_residence_time * 60),
                "{:E}".format(element.flame_depth * 0.3048),
                "{:E}".format(element.fireline_intensity_head_fire * 3.46146933),
                "{:E}".format(element.flame_length * 0.3048),
                round(element.time_of_ignition * 60, 2)
            ))


