################################################################################
##                               KP-FHR CORE                                  ##
################################################################################

from math import pi
import numpy as np
import matplotlib.pyplot as plt
import openmc
import openmc.model

################################## MATERIAL ####################################

fuel = openmc.Material(name='Fuel')
fuel.set_density('g/cm3', 10.5)
fuel.add_nuclide('U235', 4.6716e-02)
fuel.add_nuclide('U238', 2.8697e-01)
fuel.add_nuclide('O16', 5.0000e-01)
fuel.add_element('C', 1.6667e-01)

buff = openmc.Material(name='Buffer')
buff.set_density('g/cm3', 1.0)
buff.add_element('C', 1.0)
buff.add_s_alpha_beta('c_Graphite')

PyC1 = openmc.Material(name='PyC1')
PyC1.set_density('g/cm3', 1.9)
PyC1.add_element('C', 1.0)
PyC1.add_s_alpha_beta('c_Graphite')

PyC2 = openmc.Material(name='PyC2')
PyC2.set_density('g/cm3', 1.87)
PyC2.add_element('C', 1.0)
PyC2.add_s_alpha_beta('c_Graphite')

SiC = openmc.Material(name='SiC')
SiC.set_density('g/cm3', 3.2)
SiC.add_element('C', 0.5)
SiC.add_element('Si', 0.5)

graphite = openmc.Material(name='Graphite')
graphite.set_density('g/cm3', 1.1995)
graphite.add_element('C', 1.0)
graphite.add_s_alpha_beta('c_Graphite')

salt = openmc.Material(name='FLiBe')
salt.set_density('g/cm3', 1.2)
salt.add_element('F', 0.334)
salt.add_element('Li', 0.333)
salt.add_element('Be', 0.333)

structure = openmc.Material(name='SS316')
structure.set_density('g/cm3', 1.2)
structure.add_element('Fe', 0.5)
structure.add_element('C', 0.1)
structure.add_element('Cr', 0.4)

material = openmc.Materials([fuel,buff,PyC1,PyC2,SiC,graphite,salt,structure])
material.cross_sections='/home/fadilnaufal/OpenMC_CS/endfb80/cross_sections.xml'
material.export_to_xml()

################################## GEOMETRY ####################################

spheres = [openmc.Sphere(r=1e-4*d/2)
           for d in [425., 625., 705., 775.]]

cells = [openmc.Cell(fill=fuel, region=-spheres[0]),
         openmc.Cell(fill=buff, region=+spheres[0] & -spheres[1]),
         openmc.Cell(fill=PyC1, region=+spheres[1] & -spheres[2]),
         openmc.Cell(fill=SiC, region=+spheres[2] & -spheres[3]),
         openmc.Cell(fill=PyC2, region=+spheres[3])]
triso_univ = openmc.Universe(cells=cells)

fregion = [openmc.Sphere(r=df, boundary_type='transmission')
          for df in [3., 5., 6.]]

outer_radius = 855.*1e-4/2
centers=openmc.model.pack_spheres(radius=outer_radius,\
               region = +fregion[0] & -fregion[1],\
               pf=0.12)

trisos = \
[openmc.model.TRISO(outer_radius, triso_univ, center) for center in centers]
print(trisos[0])

# confirmation that all TRISO particles are within the boundaries
# centers = np.vstack([triso.center for triso in trisos])
# print(centers.min(axis=0))
# print(centers.max(axis=0))

# actual packing fraction 
# print(len(trisos)*4/3* pi * outer_radius **3)

box = openmc.Cell(region = +fregion[0] & -fregion[1])
lower_left, upper_right = box.region.bounding_box
shape = (3, 3, 3)
pitch = (upper_right - lower_left)/shape
lattice = openmc.model.create_triso_lattice(
          trisos, lower_left, pitch, shape, graphite)

box.fill = lattice

in_cell = openmc.Cell(fill=graphite, region=-fregion[0])
ou_cell = openmc.Cell(fill=graphite, region=+fregion[1] & -fregion[2])

pebb_univ = openmc.Universe(cells=[box, in_cell, ou_cell])

##########################################################

core1 = [openmc.ZCylinder(r=rc, boundary_type='transmission')
               for rc in [120., 180.]]
core2 = [openmc.ZCylinder(r=rc2, boundary_type='vacuum')
                 for rc2 in [185.]]
top = openmc.ZPlane(z0=155, boundary_type='vacuum')
bottom = openmc.ZPlane(z0=-155, boundary_type='vacuum')

core_region = -core1[0] & -top & +bottom
refl_region = +core1[0] & -core1[1] & -top & +bottom
vess_region = +core1[1] & -core2[0] & -top & +bottom

outer_radius_2 = 6.0
centers_2=openmc.model.pack_spheres(radius=outer_radius_2,\
               region = core_region,\
               pf=0.6)

pebble = \
[openmc.model.TRISO(outer_radius_2, pebb_univ, center) for center in centers_2]
print(pebble[0])

box_2 = openmc.Cell(region = core_region)
lower_left_2, upper_right_2 = box_2.region.bounding_box
shape_2 = (3, 3, 3)
pitch_2 = (upper_right_2 - lower_left_2)/shape_2
lattice_2 = openmc.model.create_triso_lattice(
          pebble, lower_left_2, pitch_2, shape_2, salt)

box_2.fill = lattice_2

refl_cell = openmc.Cell(fill=graphite, region=refl_region)
vess_cell = openmc.Cell(fill=structure, region=vess_region)

universe = openmc.Universe(cells=[box_2, refl_cell, vess_cell])

geometry = openmc.Geometry(universe)
geometry.export_to_xml()

################################## SETTINGS ####################################

settings = openmc.Settings()
settings.run_mode = 'plot'
settings.export_to_xml()

#################################### PLOT ######################################

colors = {}
colors[fuel] = 'orange'
colors[buff] = 'green'
colors[PyC1] = 'blue'
colors[SiC] = 'purple'
colors[PyC2] = 'blue'
colors[graphite] = 'gray'
colors[salt] = 'pink'
colors[structure] = 'palegreen'

plot1 = openmc.Plot.from_geometry(geometry)
plot1.basis = 'xy'
plot1.origin = (0.0, 0.0, 0.0)
plot1.width = (380, 380)
plot1.pixels = (2000, 2000)
plot1.color_by = 'material'
plot1.colors = colors
plot1.to_ipython_image()

plot2 = openmc.Plot.from_geometry(geometry)
plot2.basis = 'xz'
plot2.origin = (0.0, 0.0, 0.0)
plot2.width = (380, 380)
plot2.pixels = (2000, 2000)
plot2.color_by = 'material'
plot2.colors = colors
plot2.to_ipython_image()

p = openmc.Plots((plot1, plot2))
p.export_to_xml()
