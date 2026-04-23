```@meta
CurrentModule = ExtDocs
```

# Visualization

Plotting and figure-construction routines for reconstructed stratigraphic and facies outputs.

## High-level figures

```@docs
summary_plot
summary_plot_from_state
wheeler_diagram
sediment_profile
glamour_view
fence_plot
map_view
production_curve
stratigraphic_column
```

## Axis-mutating methods

```@docs
wheeler_diagram!
sediment_profile!
glamour_view!
production_curve!
stratigraphic_column!
stratigraphic_column_layers!
```

## Lower-level plotting helpers

```@docs
sediment_accumulation!
dominant_facies!
facies_colormap
wheeler_column
```