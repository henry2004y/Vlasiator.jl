
# Gallery

## PyPlot

* Proton density in a laminar flow with points denoting cell centers. [demo](@ref demo_plot_mesh)
![](figures/mesh.png)

* Proton density of Earth's magnetosphere in the meridional cut from 2D simulation, with fieldlines through fixed seeding points. [demo](@ref demo_2d_fieldlines)
![](figures/magnetosphere_earth_proton_density_2D_bx_bz.png)

* Proton density of Earth's magnetosphere in the meridional plane from 3D simulation.
![](figures/magnetosphere_earth_proton_density_ycut.png)

* Proton density of Earth's magnetosphere in the equatorial plane from 2D simulation, zoomed in to the magnetosheath and foreshock region, with streamlines and density contour at 10 amu/cc. [demo](@ref demo_2d_contour_streamline_levels)
![](figures/magnetosphere_earth_proton_density_2D.png)

* Proton density of Earth's magnetosphere in the normal cut planes from 3D simulation. [demo](@ref demo_3d_cuts)
![](figures/magnetosphere_earth_proton_density_3cuts.png)

* Proton phase space distribution projected onto the X-Z plane. [demo](@ref demo_vdf)
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