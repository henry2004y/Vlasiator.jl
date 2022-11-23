
# Gallery

## PyPlot

* Proton density in a laminar flow with points denoting cell centers. [demo](https://github.com/henry2004y/Vlasiator.jl/blob/master/examples/demo_mesh_pyplot.jl)
![](figures/mesh.png)

* Proton density of Earth's magnetosphere in the meridional cut from 2D simulation, with fieldlines through fixed seeding points. [demo](https://github.com/henry2004y/Vlasiator.jl/blob/master/examples/demo_B_fieldline_pyplot.jl)
![](figures/magnetosphere_earth_proton_density_2D_bx_bz.png)

* Proton density of Earth's magnetosphere in the meridional plane from 3D simulation.
![](figures/magnetosphere_earth_proton_density_ycut.png)

* Proton density of Earth's magnetosphere in the equatorial plane from 2D simulation, zoomed in to the magnetosheath and foreshock region, with streamlines and density contour at 10 amu/cc. [demo](https://github.com/henry2004y/Vlasiator.jl/blob/master/examples/demo_2dplot_pyplot.jl)
![](figures/magnetosphere_earth_proton_density_2D.png)

* Proton density of Earth's magnetosphere in the normal cut planes from 3D simulation. [demo](https://github.com/henry2004y/Vlasiator.jl/blob/master/examples/demo_3dcuts_pyplot.jl)
![](figures/magnetosphere_earth_proton_density_3cuts.png)

* Proton phase space distribution projected onto the X-Z plane. [demo](https://github.com/henry2004y/Vlasiator.jl/blob/master/examples/demo_vdf_pyplot.jl)
![](figures/phase_space_distribution.png)

## Makie

Demos can be found in the [Usage](https://github.com/henry2004y/VlasiatorMakie.jl/blob/main/README.md#usage) section of VlasiatorMakie.

* Various colored contours from 2D equatorial run
![](figures/XY_contours_makie.png)

* Interactive proton density slice from 3D AMR run
![](figures/slice_interactive.png)

* Three orthogonal slices of proton density from 3D AMR run
![](figures/3slices_makie.png)

* Isosurface of Bz = 0 from 3D AMR run
![](figures/isosurface_bz=0_makie.png)

* Proton phase space distribution projected onto the X-Z plane
![](figures/VDF_slice.png)

* Interactive proton phase space distribution in the three orthogonal planes
![](figures/VDF_slices.png)

* Proton phase space distribution
![](figures/VDF_volume.png)

## ParaView

VLSV files can be [converted to the structured VTK format](manual.md#Converting-to-VTK), and then visualized in ParaView.

* 2D slice contour of density in the meriodional plane with streamlines
![](figures/3D_paraview_slice.png)

* 2D slices of density viewing from upstream
![](figures/3D_paraview_2slices.png)