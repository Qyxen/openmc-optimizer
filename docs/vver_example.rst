Complete VVER model and optimization example
============================================

This example builds a complete VVER fuel-assembly model and uses OpenMC
Optimizer to vary uranium enrichment and Gd2O3 loading. It creates the
materials, fuel pins, hexagonal lattice, reflector, simulation settings, and
pin-power tally directly in Python; no pre-existing OpenMC XML files are
required.

The model contains one 8-ring hexagonal assembly, five fuel designs, a central
guide tube, and a reflector. Because it models one assembly rather than a full
core, a custom tally measures power in each repeated fuel-pin instance. The
optimization minimizes the pin power peaking factor (PPF).

Before running
--------------

Use a Python environment that has OpenMC, OpenMC Optimizer, and configured
OpenMC cross-section data. Run the code below as part of your Python workflow.
It creates ``materials.xml``, ``geometry.xml``, ``settings.xml``, and
``tallies.xml`` in a ``vver-model`` directory under the current working
directory. Candidate simulations and result files are stored beneath that
directory in ``vver-optimization``.
The example uses 5,000 particles per batch, 100 batches, and a population of 10
for 3 generations. Increase these settings for a production calculation.

Model and optimization source
-----------------------------

The listing below contains the complete model and optimization definition.
Its pin dimensions, lattice arrangement, materials, and starting enrichments
are specified directly. The UO2 materials use atom fractions internally, so
the code converts the requested 1–5 wt% U-235 bounds into atom fractions. The
absorber pin target materials start as pure UO2, and the separate Gd2O3 recipe
is mixed into them by two volume-fraction variables. Their bounds of 0–0.10
therefore mean 0–10 vol% Gd2O3 in each final absorber fuel.

::

   """Build the VVER model and optimize enrichment and Gd2O3 loading."""
   
   from __future__ import annotations
   
   from pathlib import Path
   
   import openmc
   
   from openmc_optimizer import MaterialFraction, OpenMCOptimizer, VolumeFraction
   
   
   MODEL_DIR = Path.cwd() / "vver-model"
   MODEL_DIR.mkdir(parents=True, exist_ok=True)
   
   
   # Uranium dioxide fuel materials at their initial enrichments
   au13 = openmc.Material(name="au13")
   au13.add_element("U", 1.0, enrichment=1.3)
   au13.add_element("O", 2.0)
   au13.set_density("g/cm3", 10.4)
   
   au22 = openmc.Material(name="au22")
   au22.add_element("U", 1.0, enrichment=2.2)
   au22.add_element("O", 2.0)
   au22.set_density("g/cm3", 10.4)
   
   av5 = openmc.Material(name="av5")
   av5.add_element("U", 1.0, enrichment=3.0)
   av5.add_element("O", 2.0)
   av5.set_density("g/cm3", 10.4)
   
   awu = openmc.Material(name="awu")
   awu.add_element("U", 1.0, enrichment=3.3)
   awu.add_element("O", 2.0)
   awu.set_density("g/cm3", 10.4)
   
   go = openmc.Material(name="go")
   go.add_element("U", 1.0, enrichment=4.0)
   go.add_element("O", 2.0)
   go.set_density("g/cm3", 10.4)
   
   gd2o3 = openmc.Material(name="Gd2O3")
   gd2o3.add_element("Gd", 2.0)
   gd2o3.add_element("O", 3.0)
   gd2o3.set_density("g/cm3", 7.41)
   
   water = openmc.Material(name="water")
   water.add_nuclide("H1", 2.0)
   water.add_nuclide("O16", 1.0)
   water.set_density("g/cm3", 1.0)
   water.add_s_alpha_beta("c_H_in_H2O")
   
   zirconium = openmc.Material(name="zirconium")
   zirconium.add_element("Zr", 1.0)
   zirconium.set_density("g/cm3", 6.6)
   
   niobium = openmc.Material(name="niobium")
   niobium.add_element("Nb", 1.0)
   niobium.set_density("g/cm3", 8.57)
   
   helium = openmc.Material(name="Helium")
   helium.add_element("He", 1.0)
   helium.set_density("g/cm3", 0.178e-3)
   
   reflector_mat = openmc.Material(name="Reflector")
   reflector_mat.add_nuclide("Be9", 1.0)
   reflector_mat.add_nuclide("O16", 1.0)
   reflector_mat.set_density("g/cm3", 2.9)
   
   
   # The absorber pin materials start as pure UO2. The optimizer adds Gd2O3.
   fuel_awu = awu
   fuel_go = go
   
   alloy = openmc.Material.mix_materials(
       [niobium, zirconium], [0.01, 0.99], "wo", name="Nb-Zr alloy"
   )
   
   materials = openmc.Materials(
       [au13, au22, av5, fuel_awu, fuel_go, gd2o3, water, alloy, helium, reflector_mat]
   )
   materials.export_to_xml(path=str(MODEL_DIR / "materials.xml"))
   
   
   # Axial and radial boundaries for the assembly and reflector.
   assembly_z0 = openmc.ZPlane(surface_id=300, z0=-75.0)
   assembly_z1 = openmc.ZPlane(surface_id=301, z0=75.0)
   reflector_z0 = openmc.ZPlane(
       surface_id=303, z0=-95.0, boundary_type="vacuum"
   )
   reflector_z1 = openmc.ZPlane(
       surface_id=304, z0=95.0, boundary_type="vacuum"
   )
   assembly_prism = -openmc.model.HexagonalPrism(edge_length=18.0, orientation="y")
   reflector_prism = -openmc.model.HexagonalPrism(
       edge_length=19.0, orientation="y", boundary_type="vacuum"
   )
   
   all_water_out = openmc.Cell(cell_id=200, fill=water)
   assembly_cell = openmc.Cell()
   reflect_cell = openmc.Cell()
   top_reflect_cell = openmc.Cell()
   bot_reflect_cell = openmc.Cell()
   
   assembly_cell.region = assembly_prism & -assembly_z1 & +assembly_z0
   reflect_cell.region = ~assembly_prism & reflector_prism & -assembly_z1 & +assembly_z0
   top_reflect_cell.region = reflector_prism & +assembly_z1 & -reflector_z1
   bot_reflect_cell.region = reflector_prism & -assembly_z0 & +reflector_z0
   reflect_cell.fill = reflector_mat
   top_reflect_cell.fill = reflector_mat
   bot_reflect_cell.fill = reflector_mat
   
   water_out_universe = openmc.Universe(cells=[all_water_out])
   
   
   # Central guide and instrument tube.
   center_fuel_surface = openmc.ZCylinder(surface_id=400, r=0.3765)
   center_clad_inner = openmc.ZCylinder(surface_id=401, r=0.4)
   center_clad_outer = openmc.ZCylinder(surface_id=402, r=0.465)
   center_fuel_region = -center_fuel_surface & -assembly_z1 & +assembly_z0
   center_clad_region = (
       +center_fuel_surface & -center_clad_outer & -assembly_z1 & +assembly_z0
   )
   center_water_region = +center_clad_outer & -assembly_z1 & +assembly_z0
   center_fuel_cell = openmc.Cell(
       cell_id=400, fill=water, region=center_fuel_region
   )
   center_clad_cell = openmc.Cell(
       cell_id=401, fill=alloy, region=center_clad_region
   )
   center_water_cell = openmc.Cell(
       cell_id=402, fill=water, region=center_water_region
   )
   central_tube_universe = openmc.Universe(
       cells=[center_fuel_cell, center_clad_cell, center_water_cell]
   )
   
   
   def make_fuel_pin(
       fuel_material: openmc.Material,
       first_id: int,
       universe_name: str,
   ) -> tuple[openmc.Universe, openmc.Cell]:
       """Build one fuel-pin universe and return it with its fuel cell."""
       fuel_surface = openmc.ZCylinder(surface_id=first_id, r=0.3765)
       clad_inner = openmc.ZCylinder(surface_id=first_id + 1, r=0.4)
       clad_outer = openmc.ZCylinder(surface_id=first_id + 2, r=0.465)
   
       fuel_region = -fuel_surface & -assembly_z1 & +assembly_z0
       gap_region = +fuel_surface & -clad_inner & -assembly_z1 & +assembly_z0
       clad_region = +clad_inner & -clad_outer & -assembly_z1 & +assembly_z0
       water_region = +clad_outer & -assembly_z1 & +assembly_z0
   
       fuel_cell = openmc.Cell(
           cell_id=first_id, fill=fuel_material, region=fuel_region
       )
       gap_cell = openmc.Cell(
           cell_id=first_id + 1, fill=helium, region=gap_region
       )
       clad_cell = openmc.Cell(
           cell_id=first_id + 2, fill=alloy, region=clad_region
       )
       water_cell = openmc.Cell(
           cell_id=first_id + 3, fill=water, region=water_region
       )
       universe = openmc.Universe(
           cells=[fuel_cell, gap_cell, clad_cell, water_cell], name=universe_name
       )
       return universe, fuel_cell
   
   
   au13_universe, au13_fuel_cell = make_fuel_pin(au13, 410, "au13 pin")
   au22_universe, au22_fuel_cell = make_fuel_pin(au22, 420, "au22 pin")
   av5_universe, av5_fuel_cell = make_fuel_pin(av5, 430, "av5 pin")
   awu_universe, awu_fuel_cell = make_fuel_pin(fuel_awu, 440, "fuel_awu pin")
   go_universe, go_fuel_cell = make_fuel_pin(fuel_go, 450, "fuel_go pin")
   
   
   # Assemble the VVER hexagonal lattice, from the outside ring inward.
   assembly_lattice = openmc.HexLattice(name="VVER fuel assembly")
   assembly_lattice.center = (0.0, 0.0)
   assembly_lattice.pitch = (1.275,)
   assembly_lattice.outer = water_out_universe
   
   ring7 = [awu_universe, awu_universe, go_universe, awu_universe,
            go_universe, awu_universe, awu_universe] * 6
   ring6 = [
       au22_universe, av5_universe, au22_universe,
       av5_universe, au22_universe, av5_universe,
       au13_universe, av5_universe, au22_universe,
       av5_universe, au22_universe, av5_universe,
       au22_universe, av5_universe, au22_universe,
       av5_universe, au22_universe, av5_universe,
       au22_universe, av5_universe, au22_universe,
       av5_universe, au22_universe, av5_universe,
       au13_universe, av5_universe, au22_universe,
       av5_universe, au22_universe, av5_universe,
       au13_universe, av5_universe, au22_universe,
       av5_universe, au22_universe, av5_universe,
   ]
   ring5 = [au13_universe, au13_universe, au22_universe,
            au13_universe, au13_universe] * 6
   ring4 = [au22_universe, av5_universe, au13_universe, av5_universe] * 6
   ring3 = [au13_universe, au22_universe, au13_universe] * 6
   ring2 = [au22_universe, av5_universe] * 6
   ring1 = [au13_universe] * 6
   ring0 = [central_tube_universe]
   assembly_lattice.universes = [ring7, ring6, ring5, ring4, ring3,
                                 ring2, ring1, ring0]
   assembly_lattice.orientation = "y"
   assembly_cell.fill = assembly_lattice
   
   root_universe = openmc.Universe(
       cells=[assembly_cell, reflect_cell, top_reflect_cell, bot_reflect_cell]
   )
   geometry = openmc.Geometry(root_universe)
   geometry.export_to_xml(path=str(MODEL_DIR / "geometry.xml"))
   
   
   # Monte Carlo simulation settings.
   settings = openmc.Settings()
   settings.batches = 100
   settings.inactive = 10
   settings.particles = 5000
   edge = 18.0
   source_space = openmc.stats.Box(
       [-edge, -edge, -edge], [edge, edge, edge]
   )
   settings.source = openmc.IndependentSource(
       space=source_space, constraints={"fissionable": True}
   )
   settings.export_to_xml(path=str(MODEL_DIR / "settings.xml"))
   
   
   # Tally power in every repeated fuel-pin instance across the assembly.
   geometry.determine_paths(instances_only=True)
   fuel_cells = [
       au13_fuel_cell,
       au22_fuel_cell,
       av5_fuel_cell,
       awu_fuel_cell,
       go_fuel_cell,
   ]
   fuel_pin_instances = []
   for fuel_cell in fuel_cells:
       fuel_pin_instances.extend(
           (fuel_cell, instance)
           for instance in range(fuel_cell.num_instances)
       )
   if len(fuel_pin_instances) < 2:
       raise ValueError("The VVER pin-power tally needs at least two fuel pins")
   
   pin_power_tally = openmc.Tally(name="vver_pin_power")
   pin_power_tally.filters = [openmc.CellInstanceFilter(fuel_pin_instances)]
   pin_power_tally.scores = ["fission-q-recoverable"]
   pin_power_tally.estimator = "tracklength"
   openmc.Tallies([pin_power_tally]).export_to_xml(
       path=str(MODEL_DIR / "tallies.xml")
   )
   
   
   # Convert weight-enrichment bounds to atom fractions for UO2 materials.
   weight_enrichment_bounds = (0.010, 0.050)
   uranium_isotopes = ("U235", "U238")
   
   
   def weight_fraction_to_atom_fraction(weight_fraction: float) -> float:
       u235_mass = 235.0439299
       u238_mass = 238.0507884
       u235_moles = weight_fraction / u235_mass
       u238_moles = (1.0 - weight_fraction) / u238_mass
       return u235_moles / (u235_moles + u238_moles)
   
   
   atom_enrichment_bounds = tuple(
       weight_fraction_to_atom_fraction(value)
       for value in weight_enrichment_bounds
   )
   
   design_variables = [
       MaterialFraction(
           material_id=au13.id,
           component="U235",
           bounds=atom_enrichment_bounds,
           name="au13_u235_enrichment",
           basis="ao",
           normalization_group=uranium_isotopes,
       ),
       MaterialFraction(
           material_id=au22.id,
           component="U235",
           bounds=atom_enrichment_bounds,
           name="au22_u235_enrichment",
           basis="ao",
           normalization_group=uranium_isotopes,
       ),
       MaterialFraction(
           material_id=av5.id,
           component="U235",
           bounds=atom_enrichment_bounds,
           name="av5_u235_enrichment",
           basis="ao",
           normalization_group=uranium_isotopes,
       ),
       MaterialFraction(
           material_id=fuel_awu.id,
           component="U235",
           bounds=atom_enrichment_bounds,
           name="fuel_awu_u235_enrichment",
           basis="ao",
           normalization_group=uranium_isotopes,
       ),
       MaterialFraction(
           material_id=fuel_go.id,
           component="U235",
           bounds=atom_enrichment_bounds,
           name="fuel_go_u235_enrichment",
           basis="ao",
           normalization_group=uranium_isotopes,
       ),
       VolumeFraction(
           target_material_id=fuel_awu.id,
           component_material_id=gd2o3.id,
           bounds=(0.0, 0.10),
           name="fuel_awu_gd2o3_volume_fraction",
       ),
       VolumeFraction(
           target_material_id=fuel_go.id,
           component_material_id=gd2o3.id,
           bounds=(0.0, 0.10),
           name="fuel_go_gd2o3_volume_fraction",
       ),
   ]
   
   optimizer = OpenMCOptimizer(
       model_dir=MODEL_DIR,
       method="genetic",
       variables=design_variables,
       ppf_tally_name="vver_pin_power",
       ppf_tally_score="fission-q-recoverable",
       target_ppf=1.1,
       target_keff=1.15,
       keff_tolerance=0.005,
       work_dir=MODEL_DIR / "vver-optimization",
       population_size=10,
       generations=3,
       elite_count=2,
       mutation_rate=0.20,
       mutation_scale=0.10,
       seed=14,
       statepoint_file="statepoint.100.h5",
   )
   
   result = optimizer.run()
   print(f"Lowest observed pin PPF: {result.best_ppf:.5f}")
   print(f"Best keff: {result.best_keff:.5f} +/- {result.best_keff_std_dev:.5f}")
   print(f"keff target met: {result.keff_target_reached}")
   print(f"PPF and keff targets met: {result.target_reached}")
   print(f"Best parameters: {result.best_parameters}")
   print(f"Summary: {result.summary_file}")
   print(f"Candidate CSV: {result.candidates_csv}")

Reading the result
------------------

The PPF is calculated as the maximum pin power divided by the mean pin power
across the tally's fuel-pin instance bins. The tally uses the
``fission-q-recoverable`` score and includes the repeated occurrences of all
five fuel-cell types. The example asks for PPF at or below 1.1 and keff within
1.15 +/- 0.005. The optimizer prefers candidates meeting the keff tolerance,
then ranks feasible candidates using its weighted PPF and keff loss. If no
candidate meets the keff tolerance, it returns the candidate closest to the
keff target. The printed status fields show whether the keff tolerance and both
targets were achieved; a reported best candidate does not imply that the goals
were reached.

``target_ppf`` is a stopping threshold as well as the requested PPF goal; it
does not replace PPF minimization. ``best_parameters`` contains U-235 atom
fractions for all five fuel designs and physical Gd2O3 volume fractions for the
two absorber fuel designs. Since their target materials start as pure UO2, each
Gd2O3 volume fraction is its actual share of the final mixture. Because this
example models one assembly, its PPF describes pins in that assembly, not a
core-wide assembly PPF.

Gd2O3 volume-fraction variables
-------------------------------

The main example includes two ``VolumeFraction`` variables. Each target ID is
one of the UO2 materials assigned to absorber fuel pins (``fuel_awu`` or
``fuel_go``); the component ID is the standalone Gd2O3 recipe. The target IDs
are the materials used by the geometry, so each candidate updates the material
that is actually in the pin. The Gd2O3 source material only needs to be present
in the material collection; it does not need to fill a cell. Both recipes have
explicit ``g/cm3`` densities so the optimizer can calculate the mixture
composition and density by physical volume.

The bounds ``(0.0, 0.10)`` request between 0 and 10 vol% Gd2O3 in the final
absorber fuel. This is the VVER-specific API setup for optimizing only these
two Gd2O3 loadings. Use it in place of the ``design_variables`` and optimizer
configuration in the full model listing; the model, XML export, geometry, and
tally setup remain as defined there::

   from openmc_optimizer import OpenMCOptimizer, VolumeFraction

   gd2o3_variables = [
       VolumeFraction(
           target_material_id=fuel_awu.id,
           component_material_id=gd2o3.id,
           bounds=(0.0, 0.10),
           name="fuel_awu_gd2o3_volume_fraction",
       ),
       VolumeFraction(
           target_material_id=fuel_go.id,
           component_material_id=gd2o3.id,
           bounds=(0.0, 0.10),
           name="fuel_go_gd2o3_volume_fraction",
       ),
   ]

   gd2o3_optimizer = OpenMCOptimizer(
       model_dir=MODEL_DIR,
       method="genetic",
       variables=gd2o3_variables,
       ppf_tally_name="vver_pin_power",
       ppf_tally_score="fission-q-recoverable",
       target_ppf=1.1,
       target_keff=1.15,
       keff_tolerance=0.005,
       work_dir=MODEL_DIR / "vver-gd2o3-optimization",
       population_size=10,
       generations=3,
       elite_count=2,
       mutation_rate=0.20,
       mutation_scale=0.10,
       seed=14,
       statepoint_file="statepoint.100.h5",
   )

   result = gd2o3_optimizer.run()
   print(f"Best PPF: {result.best_ppf:.5f}")
   print(f"Best keff: {result.best_keff:.5f} +/- {result.best_keff_std_dev:.5f}")
   print(f"Targets met: {result.target_reached}")
   print(f"Best parameters: {result.best_parameters}")

The target IDs are the material IDs actually assigned to the two absorber fuel
pins. The Gd2O3 ID names the separate source recipe in ``materials.xml``; it
does not need to fill a geometry cell. For enrichment-only optimization, use
the five ``MaterialFraction`` variables from the full listing instead. To
optimize both enrichment and Gd2O3 loading together, include both sets in the
same ``variables`` list.

Bayesian optimization variant
-----------------------------

The same VVER model and design variables can use Bayesian optimization. Replace
the genetic optimizer setup in the complete listing with this configuration.
Here ``population_size`` is the number of initial space-filling simulations;
``generations`` is the maximum number of additional Bayesian proposals. The
initial samples run first, then Bayesian optimization evaluates one proposed
candidate at a time::

   bayesian_optimizer = OpenMCOptimizer(
       model_dir=MODEL_DIR,
       method="bayesian",
       variables=design_variables,
       ppf_tally_name="vver_pin_power",
       ppf_tally_score="fission-q-recoverable",
       target_ppf=1.1,
       target_keff=1.15,
       keff_tolerance=0.005,
       work_dir=MODEL_DIR / "vver-bayesian-optimization",
       population_size=80,
       generations=20,
       seed=14,
       statepoint_file="statepoint.100.h5",
   )

   result = bayesian_optimizer.run()
   print(f"Lowest observed pin PPF: {result.best_ppf:.5f}")
   print(f"Best keff: {result.best_keff:.5f} +/- {result.best_keff_std_dev:.5f}")
   print(f"keff target met: {result.keff_target_reached}")
   print(f"PPF and keff targets met: {result.target_reached}")
   print(f"Best parameters: {result.best_parameters}")
   print(f"Summary: {result.summary_file}")
   print(f"Candidate CSV: {result.candidates_csv}")

Scikit-learn is installed with the optimizer's standard dependencies. Bayesian
optimization is useful when OpenMC evaluations are expensive: after the initial
design, it uses a Gaussian-process surrogate to select each next candidate.
