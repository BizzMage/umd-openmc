import openmc
import openmc.material
import openmc.lib
from inch_converter import in2cm

"""Editable parameters for fuel bundle position"""
fb_offset = in2cm(1.53/2) # 1.53" b/w fuel rods in 2x2 bundles
ypos = in2cm(0) # Place fuel rods 2 ft above bottom of reactor

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
    cladd_or = openmc.ZCylinder(r=in2cm(1.25 / 2 + 0.028))
    cladd_ir = openmc.ZCylinder(r=in2cm(1.25)/2)
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
    bot_fitting_cell.fill = ss304

    # Boron Carbide Interior
    b4c_cell.region = -cladd_ir & +b4c_min & -b4c_max
    b4c_cell.fill = b4c

    # Top fitting
    top_fitting_cell.region = -cladd_ir & +b4c_max & -top_fitting_max
    top_fitting_cell.fill = ss304

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

    fuel_bundle_univ = openmc.Universe(name="Fuel Bundle Universe")
    fuel_bundle_univ.add_cells(rods)
    #print(fuel_bundle_univ)
    return fuel_bundle_univ

def define_cr_fuel_bundle(fuel_universe, cr_universe, cr_position, withdraw_percent, univ_name):
    # Check for proper input
    if (withdraw_percent < 0 or withdraw_percent > 100):
        raise ValueError("Please enter a value between 0 and 100%.")

    # Control rod movement behavior
    # Bottom of fuel rod to top of uzrh fuel is 24.12125"
    max_lift = ypos + in2cm(24.12125)
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

    fuel_cr_univ = openmc.Universe(name=univ_name)
    fuel_cr_univ.add_cells(rods)
    print(fuel_cr_univ)
    return fuel_cr_univ

def define_reflector_rabbit_source():
    # Dimensions of reflector
    thickness = in2cm(2.675)
    cladd_thickness = in2cm(0.09)
    height = in2cm(24)
    # Surfaces (height here is for cross-sectional shape, which is a square)
    refl_outer = openmc.model.RectangularPrism(width=thickness + cladd_thickness, height=thickness + cladd_thickness)
    refl_inner = openmc.model.RectangularPrism(width=thickness, height=thickness)

    # Planes
    # Reflectors .31" above grid plate (assumed to be at ypos)
    refl_min = openmc.ZPlane(z0=ypos + in2cm(0.31))
    refl_max = openmc.ZPlane(z0=ypos + height + in2cm(0.31))
    # Add .5" thick ss304 on top and bottom
    refl_bot_fitting = openmc.ZPlane(z0=refl_min.z0 + in2cm(0.5))
    refl_top_fitting = openmc.ZPlane(z0=refl_max.z0 - in2cm(0.5))

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
    pass


def define_infinite_water_cell(region, name):
    cell = openmc.Cell(name=name, region=region, fill=water)
    return cell

if __name__ == "__main__":
    fuel_univ = define_fuel_rod()
    cr_univ = define_control_rod()
    fuelbundle = define_fuel_bundle(fuel_univ)
    tr_cr_univ = define_cr_fuel_bundle(fuel_univ,cr_univ,'top right',0,"Top Right CR Bundle")
    refl_univ, rabb_univ, source_univ = define_reflector_rabbit_source()

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