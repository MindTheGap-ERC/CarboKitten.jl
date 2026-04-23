```@meta
CurrentModule = ExtDocs
```

# Production and compaction

Production, compaction, preserved-layer management, and related environmental field helpers.

## Production

```@docs
production
depth_multiplier
time_multiplier
production_rate
capped_production
uniform_production
```

## Compaction primitives

```@docs
CompactionState
porosity_at_depth
compacted_thickness
compact_column!
update_sediment_height_from_compaction!
merge_layers!
```

## Layer extraction and projection

```@docs
extract_compacted_layer_maps!
project_compacted_layers_to_cube!
bury_deposition!
```


## Topography and water depth

```@docs
initial_topography
water_depth
```