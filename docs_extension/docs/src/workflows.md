```@meta
CurrentModule = ExtDocs
```

# Workflows

## From model state to reconstructed data

- [`header_from_input_exact`](@ref)
- [`datacolumn_from_state`](@ref)
- [`dataslice_from_state_exact`](@ref)
- [`datavolume_from_state`](@ref)

## From preserved layers to regular voxel grids

- [`build_vertical_coordinates`](@ref)
- [`build_selected_surfaces`](@ref)
- [`resample_to_regular_grid`](@ref)

## From reconstructed data to figures

- [`summary_plot`](@ref)
- [`summary_plot_from_state`](@ref)
- [`sediment_profile`](@ref)
- [`wheeler_diagram`](@ref)
- [`glamour_view`](@ref)
- [`fence_plot`](@ref)
- [`map_view`](@ref)
- [`production_curve`](@ref)
- [`stratigraphic_column`](@ref)

## From facies grids to interpreted environments

- [`facies_props_depthbins`](@ref)
- [`facies_props_depthbins_thickness`](@ref)
- [`facies_percent_from_layers`](@ref)
- [`facies_percent_from_ids`](@ref)
- [`classify_block`](@ref)
- [`env_selector`](@ref)
- [`build_environment_grid`](@ref)
- [`encode_environments`](@ref)
- [`write_environment_block!`](@ref)

## From simulation outputs to batch summaries

- [`run_once`](@ref)
- [`run_single_to_csv`](@ref)
- [`patch_stats_3D_from_ids`](@ref)
- [`thickness_from_topk`](@ref)
