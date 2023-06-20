from dataclasses import dataclass


@dataclass
class Fuel:
    """
    w_o:	    Oven-dry fuel load (lb/ft^2)
    ρ_p:	    Fuel particle density (lb/ft^3)
    S_τ:	    Total mineral content (fraction)
    M_x:	    Moisture of extinction (fraction)
    S_e:	    Effective mineral content (fraction)
    h:          Heat content (Btu/lb) # kj/kg to Btu/lb
    σ sigma:	Surface-area-to-volume-ratio (sav_ratio) (ft^-1)
    δ delta:	Fuel bed depth (ft)

    Fuel type:    ((dead), (live))
    Size classes: ((1h, 10h, 100h), (herb, woody))
    Sav_ratio:    ((sav_ratio_1_h, sav_ratio_10_h, sav_ratio_100_h), (sav_ratio_live_herb, sav_ratio_live_woody))
    Load:         ((load_1_h, load_10_h, load_100_h), (load_live_herb, load_live_wood))
    """
    name: str
    load: tuple
    sav_ratio: tuple
    heat_content: tuple
    bed_depth: float
    moisture_extinction_dead: float
    particle_density: float
    total_mineral_content: float
    effective_mineral_content: float
    moisture_content: tuple
    flame_length: float
