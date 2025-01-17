import openmc
import openmc.material
import openmc.lib
from inch_converter import in2cm

"""Editable parameters for fuel bundle position"""
fb_offset = in2cm(1.53/2) # 1.53" b/w fuel rods in 2x2 bundles
ypos = in2cm(19) # Place fuel rods 2 ft above bottom of reactor
facilities_offset = in2cm(11.31125) # Centerline of fuel element minus 1.64"

"""Declare material elemental makeup and cross sections."""
# Water (includes thermal scattering)
water = openmc.Material(name='Water')
water.set_density('g/cm3', 1.)
water.add_element('H', 2.)
water.add_element('O', 1.)
water.add_s_alpha_beta('c_H_in_H2O') # Thermal scattering for hydrogen in water

# 19.75% enriched uranium zirconium hydride fuel
# NOTE: wo indicates weight percent
# NOTE: cross sections include thermal scattering for H, Zr in ZrH
uzrh = openmc.Material(name='UZrH')
uzrh.set_density('g/cm3', 5.95)
uzrh.add_nuclide('U235', .01679, 'wo')
uzrh.add_nuclide('U238', .06821, 'wo')
uzrh.add_element('H', .016, 'wo')
uzrh.add_element('Zr', .899, 'wo')
""" Alt definition for enriched uranium, unused
uzrh.add_element('U', 1., enrichment=19.75)
"""

# Molybdenum
mo = openmc.Material(name='Molybdenum')
mo.set_density('g/cm3', 10.22)
mo.add_element('Mo', 1.)

# Graphite (includes thermal scattering)
graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.70)
graphite.add_element('C', 1.)
graphite.add_s_alpha_beta('c_Graphite')

# Stainless steel 
# NOTE: Numbers from OpenMC TRIGA example
ss304 = openmc.Material(name='Stainless Steel 304')
ss304.set_density('g/cm3', 8.0)
ss304.add_element('C',.002,'wo')
ss304.add_element('Si',.004,'wo')
ss304.add_element('P',.0003,'wo')
ss304.add_element('S',.0002,'wo')
ss304.add_element('V',.003,'wo')
ss304.add_element('Cr',.115,'wo')
ss304.add_element('Mn',.006,'wo')
ss304.add_element('Fe',.8495,'wo')
ss304.add_element('Ni',.005,'wo')
ss304.add_element('Mo',.01,'wo')
ss304.add_element('W',.005,'wo')

# Zirconium
zr = openmc.Material(name='Zirconium')
zr.set_density('g/cm3', 6.506)
zr.add_element('Zr', 1.)

# Air (fills empty space in fuel cells)
# NOTE: Numbers from OpenMC TRIGA example
void = openmc.Material(name='Void')
void.set_density('g/cm3', 0.001205)
void.add_element('N', 0.755268, 'wo')
void.add_element('C', 0.000124, 'wo')
void.add_element('O', 0.231781, 'wo')
void.add_element('Ar', 0.012827, 'wo')

# Aluminum
al = openmc.Material(name='Aluminum')
al.set_density('g/cm3', 2.7)
al.add_element('Al', 1.)

# Boron Carbide
b4c = openmc.Material(name='Boron Carbide')
b4c.set_density('g/cm3', 2.52)
b4c.add_element('B', 4.)
b4c.add_element('C', 1.)

# Record materials to save to model
materials = openmc.Materials([water, uzrh, mo, graphite, ss304, zr, void, al, b4c])
materials.export_to_xml()

""" Shortening code idea
# Define a helper function to parse region strings
def create_region(region_expr, plane_key, planes):
    return eval(region_expr) & +planes[f"{plane_key}_min"] & -planes[f"{plane_key}_max"]

# List of cell definitions with their regions and materials
cell_definitions = [
    (cladd_cell, "-cladd_or & +cladd_ir", "cladd", ss304),
    (bot_fitting_cell, "-bot_fitting_or", "bot_fitting", ss304),
    (bot_refl_cell, "-refl_or", "bot_refl", graphite),
    (mo_disk_cell, "-mo_disk_or", "mo_disk", mo),
    (pellet_cell, "-pellet_or & +pellet_ir", "pellet", uzrh),
    (zr_rod_cell, "-zr_rod_or", "pellet", zr),
    (top_refl_cell, "-refl_or", "top_refl", graphite),
    (air_gap_cell, "-cladd_ir", "air_gap", void),
    (top_fitting_cell, "-top_fitting_or", "top_fitting", ss304),
]

# Assign region and fill for each cell
for cell, region_expr, plane_key, material in cell_definitions:
    cell.region = create_region(region_expr, plane_key, planes)
    cell.fill = material

# Handle miscellaneous gaps between components
misc_gaps_cell.region = (
    (-cladd_ir & +refl_or & +planes["bot_refl_min"] & -planes["bot_refl_max"]) |  # Bottom reflector air gap
    (-cladd_ir & +mo_disk_or & +planes["mo_disk_min"] & -planes["mo_disk_max"]) |  # Molybdenum disk air gap
    (-pellet_ir & +zr_rod_or & +planes["pellet_min"] & -planes["pellet_max"]) |    # Air gap between Zr rod and pellets
    (-cladd_ir & +refl_or & +planes["top_refl_min"] & -planes["top_refl_max"])     # Top reflector air gap
)
misc_gaps_cell.fill = void
"""

# Assemble fuel rod
def define_fuel_rod():
    # Assert surfaces
    cladd_or = openmc.ZCylinder(r=in2cm(1.414)/2)
    top_fitting_or = openmc.ZCylinder(r=in2cm(1.37)/2) # Estimation
    bot_fitting_or = openmc.ZCylinder(r=in2cm(1.37)/2) # Estimation
    cladd_ir = openmc.ZCylinder(r=in2cm(1.37)/2)
    pellet_or = openmc.ZCylinder(r=in2cm(1.37)/2)
    mo_disk_or = openmc.ZCylinder(r=in2cm(1.366)/2)
    refl_or = openmc.ZCylinder(r=in2cm(1.291)/2)
    pellet_ir = openmc.ZCylinder(r=in2cm(0.25)/2)
    zr_rod_or = openmc.ZCylinder(r=in2cm(0.225)/2)

    # Declare lengths of each element
    bot_fitting_len = in2cm(2)
    refl_len = in2cm(3.42)
    mo_disk_len = in2cm(0.03125)
    pellet_len = in2cm(5.0*3) # 3 pellets in each fuel rod
    air_gap_len = in2cm(0.25) # Inc. from min of 0.18 to meet ref dimension
    top_fitting_len = in2cm(1.75)

    # NOTE: Bot and top fitting go .5" into cladding.
    # Cladding total length = 23.125

    # Build model from bottom up
    components = [
        ("bot_fitting", bot_fitting_len),
        ("bot_refl", refl_len),
        ("mo_disk", mo_disk_len),
        ("pellet", pellet_len),
        ("top_refl", refl_len),
        ("air_gap", air_gap_len),
        ("top_fitting", top_fitting_len),
    ]
    # Dict for component top and bottom surfaces
    planes = {}

    curr_height = 0
    for name, length in components:
        # Record plane min and max
        planes[f"{name}_min"] = openmc.ZPlane(z0=curr_height)
        curr_height = round(curr_height + length, 8)
        planes[f"{name}_max"] = openmc.ZPlane(z0=curr_height)

    # Calculate cladding min and max (0.5" inserted)
    planes["cladd_min"] = openmc.ZPlane(z0=(planes["bot_fitting_max"].z0 - in2cm(0.5)))
    planes["cladd_max"] = openmc.ZPlane(z0=(planes["top_fitting_min"].z0 + in2cm(0.50375)))
    
    # Declare cell names
    cladd_cell = openmc.Cell(name="Fuel Rod Stainless Steel Cladding")
    bot_fitting_cell = openmc.Cell(name="Fuel Rod Bottom Fitting Cell")
    bot_refl_cell = openmc.Cell(name="Fuel Rod Bottom Reflector Cell")
    mo_disk_cell = openmc.Cell(name="Fuel Rod Molybdenum Disk Cell")
    pellet_cell = openmc.Cell(name="Fuel Rod Fuel Pellets Cell")
    zr_rod_cell = openmc.Cell(name="Fuel Rod Zirconium Rod Cell")
    top_refl_cell = openmc.Cell(name="Fuel Rod Top Reflector Cell")
    air_gap_cell = openmc.Cell(name="Fuel Rod Air Gap Cell")
    top_fitting_cell = openmc.Cell(name="Fuel Rod Top Fitting Cell")
    misc_gaps_cell = openmc.Cell(name="Air fills between Fuel Cell Components Cell")
    # NOTE: Infinite water cell generated below

    # Cell regions and materials
    """Potential more concise solution:
    # List of cell definitions with their regions and materials
    cell_definitions = [
        (cladd_cell, "-cladd_or & +cladd_ir", "cladd", ss304),
        (bot_fitting_cell, "-bot_fitting_or", "bot_fitting", ss304),
        (bot_refl_cell, "-refl_or", "bot_refl", graphite),
        (mo_disk_cell, "-mo_disk_or", "mo_disk", mo),
        (pellet_cell, "-pellet_or & +pellet_ir", "pellet", uzrh),
        (zr_rod_cell, "-zr_rod_or", "pellet", zr),
        (top_refl_cell, "-refl_or", "top_refl", graphite),
        (air_gap_cell, "-cladd_ir", "air_gap", void),
        (top_fitting_cell, "-top_fitting_or", "top_fitting", ss304),
    ]

    # Assign region and fill for each cell
    for cell, region_expr, plane_key, material in cell_definitions:
        cell.region = eval(region_expr) & +planes[f"{plane_key}_min"] & -planes[f"{plane_key}_max"]
        cell.fill = material
    """
    cladd_cell.region = -cladd_or & +cladd_ir & +planes["cladd_min"] & -planes["cladd_max"]
    cladd_cell.fill = ss304

    bot_fitting_cell.region = -bot_fitting_or & +planes["bot_fitting_min"] & -planes["bot_fitting_max"]
    bot_fitting_cell.fill = ss304

    bot_refl_cell.region = -refl_or & +planes["bot_refl_min"] & -planes["bot_refl_max"]
    bot_refl_cell.fill = graphite

    mo_disk_cell.region = -mo_disk_or & +planes["mo_disk_min"] & -planes["mo_disk_max"]
    mo_disk_cell.fill = mo

    # Fuel pellet and zr rod done together (zr rod inside fuel pellets)
    pellet_cell.region = -pellet_or & +pellet_ir & +planes["pellet_min"] & -planes["pellet_max"]
    zr_rod_cell.region = -zr_rod_or & +planes["pellet_min"] & -planes["pellet_max"]
    pellet_cell.fill = uzrh
    zr_rod_cell.fill = zr

    top_refl_cell.region = -refl_or & +planes["top_refl_min"] & -planes["top_refl_max"]
    top_refl_cell.fill = graphite

    air_gap_cell.region = -cladd_ir & +planes["air_gap_min"] & -planes["air_gap_max"]
    air_gap_cell.fill = void

    top_fitting_cell.region = -top_fitting_or & +planes["top_fitting_min"] & -planes["top_fitting_max"]
    top_fitting_cell.fill = ss304

    # Create cell to handle small gaps between components
    misc_gaps_cell.region = (-cladd_ir & +refl_or & +planes["bot_refl_min"] & -planes["bot_refl_max"] | # Bottom reflector air gap b/w cladding
                            -cladd_ir & +mo_disk_or & +planes["mo_disk_min"] & -planes["mo_disk_max"] | # Molybdenum disk air gap b/w cladding
                            -pellet_ir & +zr_rod_or & +planes["pellet_min"] & -planes["pellet_max"] | # air gap b/w zr rod and uzrh pellets
                            -cladd_ir & +refl_or & +planes["top_refl_min"] & -planes["top_refl_max"]) # Top reflector air gap b/w cladding
    misc_gaps_cell.fill = void

    # Create infinite water cell surrounding the universe
    # NOTE: need to check if this is necessary
    inf_water_region = (-planes["bot_fitting_min"] | +cladd_or | +planes["top_fitting_max"] | # top, bottom, and outside cladding
                        +bot_fitting_or & +planes["bot_fitting_min"] & -planes["bot_fitting_max"] | # Bot fitting gap fill
                        +top_fitting_or & +planes["top_fitting_min"] & -planes["top_fitting_max"]) # Top fitting gap fill
    inf_water_cell = define_infinite_water_cell(inf_water_region, "Fuel Cell Infinite Surrounding Water Cell")

    fuel_univ = openmc.Universe(name="Fuel Rod Universe")
    fuel_univ.add_cells([cladd_cell, bot_fitting_cell, bot_refl_cell, mo_disk_cell,
                        pellet_cell, zr_rod_cell, top_refl_cell, air_gap_cell, top_fitting_cell,
                        misc_gaps_cell, inf_water_cell])
    print(fuel_univ)

    return fuel_univ

def define_control_rod():
    """Create control rod universe, cells and geometry
    1.25" diameter, 0.028" thick Al cladding, 17" total length. Boron carbide interior.
    Both fittings have .5" depth into cladding
    """
    # Surfaces
    cladd_or = openmc.ZCylinder(r=in2cm(1.25 / 2))
    cladd_ir = openmc.ZCylinder(r=in2cm(1.25/ 2 - 0.028))
    # Element lengths
    total_len = in2cm(17)
    cladd_len = in2cm(16)
    top_fitting_len = in2cm(1.5)
    b4c_len = in2cm(15)
    bot_fitting_len = in2cm(0.5)
    # Planes
    cladd_min = openmc.ZPlane(z0=0) # bot fitting is flush with cladding
    b4c_min = openmc.ZPlane(z0=bot_fitting_len)
    b4c_max = openmc.ZPlane(z0=(bot_fitting_len + b4c_len))
    cladd_max = openmc.ZPlane(z0=cladd_len)
    top_fitting_max = openmc.ZPlane(z0=total_len)

    # Create cells
    cladd_cell = openmc.Cell(name="Control Rod Aluminum Cladding")
    b4c_cell = openmc.Cell(name="Control Rod Boron Carbide Core")
    bot_fitting_cell = openmc.Cell(name="Control Rod Bottom Aluminum Fitting")
    top_fitting_cell = openmc.Cell(name="Control Rod Top Aluminum Fitting")
    # NOTE: Infinite water cell generated below

    # Aluminum cladding
    cladd_cell.region = -cladd_or & +cladd_ir & +cladd_min & -cladd_max
    cladd_cell.fill = al

    # Bottom Fitting
    bot_fitting_cell.region = -cladd_ir & +cladd_min & -b4c_min
    bot_fitting_cell.fill = al

    # Boron Carbide Interior
    b4c_cell.region = -cladd_ir & +b4c_min & -b4c_max
    b4c_cell.fill = b4c

    # Top fitting
    top_fitting_cell.region = -cladd_ir & +b4c_max & -top_fitting_max
    top_fitting_cell.fill = al

    # Infinite surrounding water
    inf_water_region = (+cladd_or | +top_fitting_max | -cladd_min | # top, bottom, and outside cladding
                       +cladd_ir & +cladd_max & -top_fitting_max) # Gap between top fitting and cladding
    inf_water_cell = define_infinite_water_cell(inf_water_region, "Control Rod Surrounding Infinite Water Cell")

    cr_univ = openmc.Universe(name="Control Rod Universe")
    cr_univ.add_cells([cladd_cell, bot_fitting_cell, b4c_cell, top_fitting_cell, inf_water_cell])
    print(cr_univ)

    return cr_univ

def define_fuel_bundle(fuel_universe):
    # Outside dimension of fuel cladding
    radius = in2cm(1.414/2)
    # Declare rod positions
    rod_positions = [
        ("Top Left Rod", (-fb_offset, fb_offset, ypos)),
        ("Top Right Rod", (fb_offset, fb_offset, ypos)),
        ("Bottom Left Rod", (-fb_offset, -fb_offset, ypos)),
        ("Bottom Right Rod", (fb_offset, -fb_offset, ypos)),
    ]

    # Loop through rods to define regions
    rods = []
    regions = []

    for name, translation in rod_positions:
        # Create cell and region
        rod = openmc.Cell(name=name, fill=fuel_universe)
        rod.translation = translation

        region = openmc.ZCylinder(x0=translation[0], y0=translation[1], r=radius)
        rod.region = -region

        rods.append(rod)
        regions.append(region)

    inf_water_region = +regions[0] & +regions[1] & +regions[2] & +regions[3]
    inf_water_cell = define_infinite_water_cell(inf_water_region, name="Fuel Bundle Infinite Water")
    rods.append(inf_water_cell)

    fuel_bundle_univ = openmc.Universe(name="Fuel Bundle Universe")
    fuel_bundle_univ.add_cells(rods)
    #print(fuel_bundle_univ)
    return fuel_bundle_univ

def define_cr_fuel_bundle(fuel_universe, cr_universe, cr_position, withdraw_percent, univ_name):
    # Check for proper input
    if (withdraw_percent < 0 or withdraw_percent > 100):
        raise ValueError("Please enter a value between 0 and 100%.")

    # Control rod movement behavior
    # Bottom of fuel rod to top of uzrh fuel is 20.45125". Remove 0.5" from bottom of CR that is just aluminum
    max_lift = in2cm(20.45125 - 0.5)
    max_drop = in2cm(15) # cr drops 15 inches to bottom position
    cr_lift = max_lift - max_drop * (1 - withdraw_percent/100)
    # Outside dimension of fuel cladding
    radius = in2cm(1.414/2)

    # Place control rod
    match cr_position.lower():
        case "top left":
            rod_positions = [
                ("Top Left Control Rod", (-fb_offset, fb_offset, ypos + cr_lift), cr_universe),
                ("Top Right Rod", (fb_offset, fb_offset, ypos), fuel_universe),
                ("Bottom Left Rod", (-fb_offset, -fb_offset, ypos), fuel_universe),
                ("Bottom Right Rod", (fb_offset, -fb_offset, ypos), fuel_universe),
            ]
        case "top right":
            rod_positions = [
                ("Top Left Rod", (-fb_offset, fb_offset, ypos), fuel_universe),
                ("Top Right Control Rod", (fb_offset, fb_offset, ypos + cr_lift), cr_universe),
                ("Bottom Left Rod", (-fb_offset, -fb_offset, ypos), fuel_universe),
                ("Bottom Right Rod", (fb_offset, -fb_offset, ypos), fuel_universe),
            ]
        case _:
            raise ValueError("Incorrect cr position. Please enter 'top right/left'")

    # Loop through rods to define regions
    rods = []
    regions = []

    for name, translation, universe in rod_positions:
        # Create cell and region
        rod = openmc.Cell(name=name, fill=universe)
        rod.translation = translation

        region = openmc.ZCylinder(x0=translation[0], y0=translation[1], r=radius)
        rod.region = -region

        rods.append(rod)
        regions.append(region)

    inf_water_region = +regions[0] & +regions[1] & +regions[2] & +regions[3]
    inf_water_cell = define_infinite_water_cell(inf_water_region, name="Control Rod Infinite Water")
    rods.append(inf_water_cell)

    fuel_cr_univ = openmc.Universe(name=univ_name)
    fuel_cr_univ.add_cells(rods)
    print(fuel_cr_univ)
    return fuel_cr_univ

def define_reflector_rabbit_source():
    # Dimensions of reflector
    thickness = in2cm(2.775)
    cladd_thickness = in2cm(0.016)
    height = in2cm(25)
    # Surfaces (height here is for cross-sectional shape, which is a square)
    refl_outer = openmc.model.RectangularPrism(width=thickness + cladd_thickness, height=thickness + cladd_thickness)
    refl_inner = openmc.model.RectangularPrism(width=thickness, height=thickness)

    # Planes
    # Reflectors .31" above grid plate (assumed to be at ypos)
    refl_min = openmc.ZPlane(z0=ypos + in2cm(0.31))
    refl_max = openmc.ZPlane(z0=ypos + height + in2cm(0.31))
    # Add .5" thick ss304 on top and bottom
    refl_bot_fitting = openmc.ZPlane(z0=refl_min.z0 + in2cm(0.065))
    refl_top_fitting = openmc.ZPlane(z0=refl_max.z0 - in2cm(0.065))

    # Create cells
    refl_cladd_cell = openmc.Cell(name="Graphite Reflector Aluminum Cladding Cell")
    refl_core_cell = openmc.Cell(name="Graphite Reflector Core Cell")
    refl_bot_fitting_cell = openmc.Cell(name="Graphite Reflector Bottom Aluminum Fitting Cell")
    refl_top_fitting_cell = openmc.Cell(name="Graphite Reflector Top Aluminum Fitting Cell")
    # Infinite water cell generated below

    refl_cladd_cell.region = -refl_outer & +refl_inner & +refl_min & -refl_max
    refl_cladd_cell.fill = al

    refl_core_cell.region = -refl_inner & +refl_bot_fitting & -refl_top_fitting
    refl_core_cell.fill = graphite

    refl_bot_fitting_cell.region = -refl_inner & +refl_min & -refl_bot_fitting
    refl_top_fitting_cell.region = -refl_inner & +refl_top_fitting & -refl_max
    refl_bot_fitting_cell.fill = al
    refl_top_fitting_cell.fill = al

    inf_water_region = +refl_max | -refl_min | +refl_outer
    inf_water_cell = define_infinite_water_cell(inf_water_region, "Reflector Infinite Surrounding Water")

    refl_univ = openmc.Universe(name="Graphite Reflector Universe")
    refl_univ.add_cells([refl_cladd_cell, refl_core_cell, refl_top_fitting_cell, refl_bot_fitting_cell, inf_water_cell])
    #print(refl_univ)

    ### Rabbit shares same geometry, just with an air hole in the center
    rabb_center_hole = openmc.ZCylinder(r=in2cm(1.75)/2)
    # Cell names
    rabb_cladd_cell = openmc.Cell(name="Rabbit Aluminum Cladding Cell")
    rabb_core_cell = openmc.Cell(name="Rabbit Graphite Core Cell")
    rabb_hole_cell = openmc.Cell(name="Rabbit Air Hole Cell")
    rabb_bot_fitting_cell = openmc.Cell(name="Rabbit Bottom Aluminum Fitting Cell")
    rabb_top_fitting_cell = openmc.Cell(name="Rabbit Top Aluminum Fitting Cell")
    # Infinite water generated further down

    rabb_cladd_cell.region = -refl_outer & +refl_inner & +refl_min & -refl_max
    rabb_cladd_cell.fill = al

    rabb_core_cell.region = -refl_inner & +rabb_center_hole & +refl_bot_fitting & -refl_top_fitting
    rabb_core_cell.fill = graphite

    rabb_hole_cell.region = -rabb_center_hole & +refl_bot_fitting & -refl_max
    rabb_hole_cell.fill = void

    rabb_bot_fitting_cell.region = -refl_inner & +refl_min & -refl_bot_fitting
    rabb_top_fitting_cell.region = -refl_inner & +rabb_center_hole & +refl_top_fitting & -refl_max
    rabb_bot_fitting_cell.fill = al
    rabb_top_fitting_cell.fill = al

    rabb_inf_water_cell = define_infinite_water_cell(inf_water_region, "Rabbit Surrounding Infinite Water Cell")

    rabb_univ = openmc.Universe(name="Rabbit Universe")
    rabb_univ.add_cells([rabb_cladd_cell, rabb_core_cell, rabb_hole_cell, rabb_bot_fitting_cell, rabb_top_fitting_cell, rabb_inf_water_cell])

    ### PuBe source holder also shares reflector geometry, so it's in here too
    source_center_hole = openmc.ZCylinder(r=in2cm(1/2)) # 1" center hole filled with water

    # Cell names
    source_cladd_cell = openmc.Cell(name="PuBe Source Aluminum Cladding Cell")
    source_core_cell = openmc.Cell(name="PuBe Source Core Cell")
    source_hole_cell = openmc.Cell(name="PuBe Source Water Hole Cell")
    source_bot_fitting_cell = openmc.Cell(name="PuBe Source Bottom Aluminum Fitting Cell")
    source_top_fitting_cell = openmc.Cell(name="PuBe Source Top Aluminum Fitting Cell")

    source_cladd_cell.region = -refl_outer & +refl_inner & +refl_min & -refl_max
    source_cladd_cell.fill = al

    source_core_cell.region = -refl_inner & +source_center_hole & +refl_bot_fitting & -refl_top_fitting
    source_core_cell.fill = graphite

    source_hole_cell.region = -source_center_hole & +refl_bot_fitting & -refl_max
    source_hole_cell.fill = water

    source_bot_fitting_cell.region = -refl_inner & +refl_min & -refl_bot_fitting
    source_top_fitting_cell.region = -refl_inner & +source_center_hole & +refl_top_fitting & -refl_max
    source_bot_fitting_cell.fill = al
    source_top_fitting_cell.fill = al

    source_inf_water_cell = define_infinite_water_cell(inf_water_region, "PuBe Source Surrounding Infinite Water Cell")

    source_univ = openmc.Universe(name="PuBe Source Universe")
    source_univ.add_cells([source_cladd_cell, source_core_cell, source_hole_cell, source_bot_fitting_cell, source_top_fitting_cell, source_inf_water_cell])
    #print(source_univ)

    ### Return all the universes
    return refl_univ, rabb_univ, source_univ

def define_thermal_column():
    width = in2cm(24)
    cladd_width = in2cm(0.25)
    #24x24" and 60" thick
    tc_outer = openmc.model.RectangularPrism(width=width+cladd_width, height=width+cladd_width, axis='y')
    tc_inner = openmc.model.RectangularPrism(width=width, height=width, axis='y')

    # planes
    tc_min = openmc.YPlane(y0=0)
    tc_max = openmc.YPlane(y0=in2cm(60.5))
    # .5" thick aluminum on front and back
    tc_front_fitting = openmc.YPlane(y0=in2cm(0.25))
    tc_back_fitting = openmc.YPlane(y0=tc_max.y0 - in2cm(0.25))

    # Cell names
    tc_cladd_cell = openmc.Cell(name="Thermal Column Cladding Cell")
    tc_core_cell = openmc.Cell(name="Thermal Column Graphite Core Cell")
    tc_front_fitting_cell = openmc.Cell(name="Thermal Column Bottom Aluminum Fitting Cell")
    tc_back_fitting_cell = openmc.Cell(name="Thermal Column Top Aluminum Fitting Cell")
    # Infinite water cell below

    tc_cladd_cell.region = -tc_outer & +tc_inner & +tc_min & -tc_max
    tc_cladd_cell.fill = al

    tc_core_cell.region = -tc_inner & +tc_front_fitting & -tc_back_fitting
    tc_core_cell.fill = graphite

    tc_front_fitting_cell.region = -tc_inner & +tc_min & -tc_front_fitting
    tc_front_fitting_cell.fill = al

    tc_back_fitting_cell.region = -tc_inner & +tc_back_fitting & -tc_max
    tc_back_fitting_cell.fill = al

    inf_water_region = +tc_outer | -tc_min | +tc_max
    inf_water_cell = define_infinite_water_cell(inf_water_region, "Thermal Column Surrounding Infinite Water Cell")

    tc_univ = openmc.Universe(name="Thermal Column Universe")
    tc_univ.add_cells([tc_cladd_cell, tc_core_cell, tc_front_fitting_cell, tc_back_fitting_cell, inf_water_cell])
    #print(tc_univ)
    return tc_univ

def define_tubes():
    # 5.875" OD, 5.375" ID air filled aluminum tube that is centered 3.25" north of the north face of the core.
    through_tube_od = openmc.XCylinder(r=in2cm(5.875) / 2) # Says diameter but this only takes in radius
    through_tube_id = openmc.XCylinder(r=in2cm(5.375) / 2)
    beam_end = openmc.XPlane(x0=0) # Used later for east and west tubes

    cell_properties = [
        {"name": "Through Tube Aluminum Cladding Cell", "region": -through_tube_od & +through_tube_id, "fill": al},
        {"name": "Through Tube Air Cell", "region": -through_tube_id, "fill": void},
        {"name": "Through Tube Infinite Surrounding Water Cell", "region": +through_tube_od, "fill": water},
    ]

    cells = [openmc.Cell(name=cell["name"], region=cell["region"], fill=cell["fill"]) for cell in cell_properties]
    
    through_tube_univ = openmc.Universe(name="Through Tube Universe")
    through_tube_univ.add_cells(cells)
    #print(through_tube_univ)

    ### Define the west tube (positive x direction)
    west_beam_endcap = openmc.XPlane(x0=in2cm(0.5))
    west_beam_region = (+beam_end & -through_tube_od & +through_tube_id | # Tube walls
                        +beam_end & -west_beam_endcap & -through_tube_id) # Tube endcap
    west_cell_properties = [
        {"name": "West Tube Aluminum Cladding Cell", "region": west_beam_region, "fill": al},
        {"name": "West Tube Air Cell", "region": +west_beam_endcap & -through_tube_id, "fill": void},
        {"name": "West Tube Infinite Surrounding Water Cell", "region": -beam_end | +through_tube_od, "fill": water},
    ]

    west_beam_cells = [openmc.Cell(name=cell["name"], region=cell["region"], fill=cell["fill"]) for cell in west_cell_properties]

    west_beam_univ = openmc.Universe(name="West Beam Port Universe")
    west_beam_univ.add_cells(west_beam_cells)

    ### Now the same for east beam cells
    east_beam_endcap = openmc.XPlane(x0=-in2cm(0.5))
    east_beam_region = (-beam_end & -through_tube_od & +through_tube_id | # Tube walls
                        -beam_end & +east_beam_endcap & -through_tube_id) # Tube endcap
    east_cell_properties = [
        {"name": "East Tube Aluminum Cladding Cell", "region": east_beam_region, "fill": al},
        {"name": "East Tube Air Cell", "region": -east_beam_endcap & -through_tube_id, "fill": void},
        {"name": "East Tube Infinite Surrounding Water Cell", "region": +beam_end | +through_tube_od, "fill": water},
    ]

    east_beam_cells = [openmc.Cell(name=cell["name"], region=cell["region"], fill=cell["fill"]) for cell in east_cell_properties]

    east_beam_univ = openmc.Universe(name="East Beam Port Universe")
    east_beam_univ.add_cells(east_beam_cells)

    return through_tube_univ, west_beam_univ, east_beam_univ

def define_lattice(f, l, r, c, e, a, p, config):
    """ Order of inputs:
        Fuel, shim1, shim2, regrod, reflector, rabbit, source
    """
    # Pitch for lattice (distance from center to center for each bundle)
    xpitch = in2cm(3.189)
    ypitch = in2cm(3.035)

    # Empty cells for water gaps in lattice
    water_pin_cell = openmc.Cell(name="Empty Water Pin Cell for Lattice", fill=water)
    lattice_inf_cell = openmc.Cell(name="Infinite Surrounding Water Cell Outside Lattice", fill=water)

    # water universe name shortened for lattice
    w = openmc.Universe(name="Empty Water Pin Universe for Lattice", cells=(water_pin_cell,))
    lattice_inf_univ = openmc.Universe(name="Infinite Surrounding Water Universe Outside Lattice", cells=(lattice_inf_cell,))

    # Lattice properties
    lower_left = [-4.5*xpitch, -2.5*ypitch]
    pitch = [xpitch, ypitch]
    outer_univ = lattice_inf_univ

    lattice = openmc.RectLattice(name="Fuel Rod Lattice Structure")
    lattice.lower_left = lower_left
    lattice.pitch = pitch
    lattice.outer = outer_univ

    match config.lower():
        case "modern":
            lattice.universes = [[w, f, f, f, f, f, f, w, w],
                                [w, f, r, f, f, l, f, e, w],
                                [w, f, f, f, f, f, f, e, w],
                                [w, f, f, c, f, a, f, w, w],
                                [w, f, f, p, f, f, f, w, w]]   
        case "old":
            lattice.universes = [[w, f, f, f, f, f, f, w, w],
                                [w, f, r, f, f, l, f, w, w],
                                [w, f, f, f, f, f, f, w, w],
                                [w, f, f, c, f, a, f, w, w],
                                [w, w, w, p, e, f, e, w, w]]
        case _:
            raise ValueError("Improper reactor config. Enter modern/old")

    return lattice

def define_infinite_water_cell(region, name):
    cell = openmc.Cell(name=name, region=region, fill=water)
    return cell

def assemble_root_universe(lattice, through_tube, west_beam, east_beam, thermal_col):
    # Tube radius
    tube_radius = in2cm(5.875) / 2
    # xpitch, ypitch from above
    xpitch = in2cm(3.189)
    ypitch = in2cm(3.035)
    n_s_core_edge = ypitch*2.5
    e_w_core_edge = xpitch*4.5
    # Place lattice
    lattice_cell_area = openmc.model.RectangularPrism(width=xpitch*9, height=ypitch*5, axis='z')

    lattice_cell = openmc.Cell(name="Lattice Cell", region=-lattice_cell_area, fill=lattice)

    # Through tube
    # centered 3.25" north of the north face of the core
    through_tube_area = openmc.XCylinder(y0= -n_s_core_edge -in2cm(3.25), z0=ypos+facilities_offset, r=tube_radius)
    through_tube_cell = openmc.Cell(name="Through Tube Cell", region = -through_tube_area, fill=through_tube)
    through_tube_cell.translation = (0, through_tube_area.y0, through_tube_area.z0)

    # West Beam Port
    # ends 4.5" from the west side of the reactor core (including the reflectors). 
    # It is centered 7.7" north of south face of the core
    west_beam_area = openmc.XCylinder(y0=n_s_core_edge - in2cm(7.7), z0=ypos+facilities_offset, r=tube_radius)
    west_beam_min = openmc.XPlane(x0=e_w_core_edge + in2cm(4.5))
    west_beam_cell = openmc.Cell(name="West Beam Port Cell", region= +west_beam_min & -west_beam_area, fill= west_beam)
    west_beam_cell.translation = (west_beam_min.x0, west_beam_area.y0, west_beam_area.z0)

    # East Beam Port
    # ends 3.8" from the east side of the reactor core. 
    # It is centered 7.7" north of south face of the core. 
    east_beam_area = openmc.XCylinder(y0=west_beam_area.y0, z0=west_beam_area.z0, r=tube_radius)
    east_beam_max = openmc.XPlane(x0= -e_w_core_edge - in2cm(3.8))
    east_beam_cell = openmc.Cell(name="East Beam Port Cell", region = -east_beam_max & -east_beam_area, fill=east_beam)
    east_beam_cell.translation = (east_beam_max.x0, east_beam_area.y0, east_beam_area.z0)

    # Thermal Column
    # 24x24x60 located 1" from south of core, .25" thick cladding
    thermal_col_area = openmc.model.RectangularPrism(width=in2cm(24.5), height=in2cm(24.5), axis='y', origin=(0,ypos+facilities_offset))
    thermal_col_min = openmc.YPlane(y0=n_s_core_edge + in2cm(1))
    thermal_col_max = openmc.YPlane(y0=thermal_col_min.y0 + in2cm(60.5))
    thermal_col_cell = openmc.Cell(name="Thermal Column Cell", region=-thermal_col_area & +thermal_col_min & -thermal_col_max, fill=thermal_col)
    thermal_col_cell.translation = (0, thermal_col_min.y0, ypos+facilities_offset)

    # Infinnite Surrounding Water
    inf_reactor_cell = openmc.Cell(name="Infinite Water Reactor Pool Cell", fill=water)

    # Create universe
    inf_reactor_univ = openmc.Universe(name="infinite Reactor Pool Universe")
    inf_reactor_univ.add_cells([lattice_cell, through_tube_cell,west_beam_cell, east_beam_cell, thermal_col_cell, inf_reactor_cell])
    
    # Place infinite reactor into finite reactor
    # 7ft diameter, 20'6" deep
    pool_ir = openmc.ZCylinder(r=in2cm(7*12), boundary_type='vacuum')
    pool_min = openmc.ZPlane(z0=0, boundary_type='vacuum')
    pool_max = openmc.ZPlane(z0=in2cm(20*12 + 6), boundary_type='vacuum')

    interior = openmc.Cell(name="Reactor Interior", region= -pool_ir & +pool_min & -pool_max,fill=inf_reactor_univ)

    reactor_root_univ = openmc.Universe(name="Reactor Root Universe", cells=(interior,))
    #print(reactor_root_univ)
    return reactor_root_univ

import openmc

import openmc

def plot_universe(universe, axes="xy", offset=(0.0, 0.0, 0.0)):
    """
    Plots a 2D cross-section of a given universe with auto-centered origin and offset.
    
    Parameters:
        universe (openmc.Universe): The universe to plot.
        axes (str): The plane to plot ("xy", "xz", or "yz").
        offset (tuple): Offset to apply to the center (dx, dy, dz).
    """
    # Export geometry of universe
    geometry = openmc.Geometry(universe)
    geometry.export_to_xml()
    # Get the bounding box of the universe
    bounding_box = universe.bounding_box
    lower_left = bounding_box[0]
    upper_right = bounding_box[1]

    # Calculate the center of the bounding box
    center = [
        (lower_left[i] + upper_right[i]) / 2 for i in range(3)
    ]
    
    # Apply the offset
    origin = [center[i] + offset[i] for i in range(3)]

    # Determine width of the plot based on bounding box size
    width = [
        abs(upper_right[i] - lower_left[i]) for i in range(3)
    ]

    # Set plot dimensions and basis
    if axes == "xy":
        plot_width = (width[0], width[1])
        plot_origin = (origin[0], origin[1], origin[2])  # z is for depth
        basis = "xy"
    elif axes == "xz":
        plot_width = (width[0], width[2])
        plot_origin = (origin[0], origin[1], origin[2])  # y is for depth
        basis = "xz"
    elif axes == "yz":
        plot_width = (width[1], width[2])
        plot_origin = (origin[0], origin[1], origin[2])  # x is for depth
        basis = "yz"
    else:
        raise ValueError(f"Invalid axes '{axes}'. Choose from 'xy', 'xz', 'yz'.")

    # Create the OpenMC plot object
    plot = openmc.Plot()
    plot.width = plot_width
    plot.origin = plot_origin
    plot.basis = basis
    plot.color_by = "material"  # Can also use 'cell' if desired

    # Create the Plots collection and export to XML
    plots = openmc.Plots([plot])
    plots.export_to_xml()

    # Optionally, generate the plot
    openmc.plot_geometry()

if __name__ == "__main__":
    shim1_withdraw = 90.6
    shim2_withdraw = 94.1
    regrod_withdraw = 0

    fuel_univ = define_fuel_rod()
    cr_univ = define_control_rod()
    fuelbundle = define_fuel_bundle(fuel_univ)
    shim1 = define_cr_fuel_bundle(fuel_univ, cr_univ, 'top left', shim1_withdraw, "Shim 1")
    shim2 = define_cr_fuel_bundle(fuel_univ, cr_univ, 'top right', shim2_withdraw, "Shim 2")
    regrod = define_cr_fuel_bundle(fuel_univ, cr_univ, 'top right', regrod_withdraw, "Reg Rod")
    refl_univ, rabb_univ, source_univ = define_reflector_rabbit_source()

    tc_univ = define_thermal_column()
    tube_univ, west_beam_univ, east_beam_univ = define_tubes()
    lattice = define_lattice(fuelbundle, shim1, shim2, regrod, refl_univ, rabb_univ, source_univ, 'modern')
    root_univ = assemble_root_universe(lattice, tube_univ, west_beam_univ, east_beam_univ, tc_univ)

    #plot_universe(root_univ, 'xy')
    geometry = openmc.Geometry(regrod)
    geometry.export_to_xml()
    cr_plot = openmc.Plot()
    cr_plot.basis = 'xz'
    cr_plot.origin = [0, fb_offset, 90]
    cr_plot.color_by = 'material'
    cr_plot.pixels = (4000,4000)
    cr_plot.width=(80,80)
    
    plots = openmc.Plots([cr_plot])
    plots.export_to_xml()
    openmc.plot_geometry()

    geometry = openmc.Geometry(shim1)
    geometry.export_to_xml()
    cr_plot = openmc.Plot()
    cr_plot.basis = 'xz'
    cr_plot.origin = [0, fb_offset, 90]
    cr_plot.color_by = 'material'
    cr_plot.pixels = (4000,4000)
    cr_plot.width=(80,80)
    
    plots = openmc.Plots([cr_plot])
    plots.export_to_xml()
    openmc.plot_geometry()

    geometry = openmc.Geometry(shim2)
    geometry.export_to_xml()
    cr_plot = openmc.Plot()
    cr_plot.basis = 'xz'
    cr_plot.origin = [0, fb_offset, 90]
    cr_plot.color_by = 'material'
    cr_plot.pixels = (4000,4000)
    cr_plot.width=(80,80)
    
    plots = openmc.Plots([cr_plot])
    plots.export_to_xml()
    openmc.plot_geometry()

    geometry = openmc.Geometry(root_univ)
    geometry.export_to_xml()

    cr_plot = openmc.Plot()
    cr_plot.basis = 'xz'
    cr_plot.origin = [0, -6, 90]
    cr_plot.color_by = 'material'
    cr_plot.pixels = (4000,4000)
    cr_plot.width=(80,80)
    
    plots = openmc.Plots([cr_plot])
    plots.export_to_xml()
    openmc.plot_geometry()

    settings = openmc.Settings()
    settings.particles = 50000
    settings.batches = 200
    settings.inactive = 50
    settings.run_mode = 'eigenvalue'


    # Define source
    bounds = [-36.45027, -19.27225, 48.26, 36.45027, 19.27225, 124.46]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source = openmc.IndependentSource(
        space=uniform_dist,
        constraints={'fissionable': True}
        )

    geometry = openmc.Geometry(root_univ)
    geometry.export_to_xml()
    
    settings.export_to_xml()
    #openmc.run()

    """
    print(root_univ)

    geometry = openmc.Geometry(shim1)
    geometry.export_to_xml()
    cr_plot = openmc.Plot()
    cr_plot.basis = 'yz'
    cr_plot.origin = [-fb_offset, 0, 90]
    cr_plot.color_by = 'material'
    cr_plot.pixels = (400,400)
    cr_plot.width=(80,80)
    
    plots = openmc.Plots([cr_plot])
    plots.export_to_xml()
    openmc.plot_geometry()

    geometry = openmc.Geometry(root_univ)
    geometry.export_to_xml()
    fb_plot = openmc.Plot()
    fb_plot.basis = 'xy'
    fb_plot.origin = [0, 0, 76.990575]
    fb_plot.color_by = 'material'
    fb_plot.width=(100,100)
    fb_plot.pixels=(5000,5000)
    plots = openmc.Plots([fb_plot])
    plots.export_to_xml()
    openmc.plot_geometry()
    """

  

    """
    fb_plot = openmc.Plot()
    fb_plot.basis = 'xy'
    fb_plot.origin = [0, 0, 30]
    fb_plot.color_by = 'material'
    fb_plot.width=(100,100)
    fb_plot.pixels=(5000,5000)
    plots = openmc.Plots([fb_plot])
    plots.export_to_xml()
    openmc.plot_geometry()

    geometry = openmc.Geometry(cr_univ)
    geometry.export_to_xml()

    cr_plot = openmc.Plot()
    cr_plot.basis = 'yz'
    cr_plot.origin = [0, 0, 22.5]
    cr_plot.color_by = 'material'
    cr_plot.pixels = (400,400)
    cr_plot.width=(10,80)
    
    plots = openmc.Plots([cr_plot])
    plots.export_to_xml()
    openmc.plot_geometry()
    """

    """
    geometry = openmc.Geometry(refl_univ)
    geometry.export_to_xml()
    fuel_plot = openmc.Plot()
    fuel_plot.basis = 'yz'
    fuel_plot.origin = [0, 0, 32.5]
    fuel_plot.color_by = 'material'
    fuel_plot.colors = {void:'cyan'}
    fuel_plot.width=(10,70)

    plots = openmc.Plots([fuel_plot])
    plots.export_to_xml()
    openmc.plot_geometry()

    geometry = openmc.Geometry(rabb_univ)
    geometry.export_to_xml()
    fuel_plot = openmc.Plot()
    fuel_plot.basis = 'yz'
    fuel_plot.origin = [0, 0, 32.5]
    fuel_plot.color_by = 'material'
    fuel_plot.colors = {void:'cyan'}
    fuel_plot.width=(10,70)

    plots = openmc.Plots([fuel_plot])
    plots.export_to_xml()
    openmc.plot_geometry()

    geometry = openmc.Geometry(source_univ)
    geometry.export_to_xml()
    fuel_plot = openmc.Plot()
    fuel_plot.basis = 'yz'
    fuel_plot.origin = [0, 0, 32.5]
    fuel_plot.color_by = 'material'
    fuel_plot.colors = {void:'cyan'}
    fuel_plot.width=(10,70)

    plots = openmc.Plots([fuel_plot])
    plots.export_to_xml()
    openmc.plot_geometry()
    """


    """
    fb_plot = openmc.Plot()
    fb_plot.basis = 'xy'
    fb_plot.origin = [0, 0, 32.5]
    fb_plot.color_by = 'material'
    fb_plot.colors = {b4c:'cyan'}
    fb_plot.width=(10,10)
    


    CR Plots
    geometry = openmc.Geometry(cr_univ)
    geometry.export_to_xml()

    cr_plot = openmc.Plot()
    cr_plot.basis = 'yz'
    cr_plot.origin = [0, 0, 22.5]
    cr_plot.color_by = 'material'
    cr_plot.width=(10,50)
    
    plots = openmc.Plots([cr_plot])
    plots.export_to_xml()
    openmc.plot_geometry()
    """
    # Break
    """ Fuel Universe Plots
    geometry = openmc.Geometry(fuel_univ)
    geometry.export_to_xml()
    fuel_plot = openmc.Plot()
    fuel_plot.basis = 'yz'
    fuel_plot.origin = [0, 0, 32.5]
    fuel_plot.color_by = 'material'
    fuel_plot.colors = {void:'cyan'}
    fuel_plot.width=(10,70)

    mo_plot = openmc.Plot()
    mo_plot.basis = 'yz'
    mo_plot.origin = [0, 0, 14]
    mo_plot.color_by = 'material'
    mo_plot.colors = {void:'cyan'}
    mo_plot.width=(4,10)

    plots = openmc.Plots([fuel_plot, mo_plot])
    plots.export_to_xml()
    openmc.plot_geometry()
    """