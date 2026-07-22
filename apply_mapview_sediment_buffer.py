#!/usr/bin/env python3
"""Patch CarboKitten PR #244 to reconstruct map-view layers with SedimentStack.

Run this script from the CarboKitten.jl repository root on branch
243-Add-map-view-visualization.

The patch is deliberately additive and backward compatible:
- existing `map_view` / `map_view!` calls remain valid;
- `show = :model` is unchanged;
- `show = :preserved` and `show = :both` reconstruct all saved deposition and
  disintegration up to the selected time with the existing sediment stack;
- optional `depth`, `layer_thickness`, and `depositional_resolution` keywords
  control the sampled preserved interval;
- only the Entangled Markdown sources are edited directly. Generated Julia
  files are refreshed with `entangled tangle`, avoiding a stitch/tangle conflict.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path
from typing import Callable, NoReturn

BRANCH = "243-Add-map-view-visualization"

SEDIMENT_DOC = Path("docs/src/components/sediment_buffer.md")
SEDIMENT_SRC = Path("src/SedimentStack.jl")
SEDIMENT_TEST = Path("test/SedimentStackSpec.jl")
MAP_DOC = Path("docs/src/visualization/map-view.md")
MAP_SRC = Path("ext/MapView.jl")


def fail(message: str) -> NoReturn:
    raise SystemExit(message)


def current_branch() -> str:
    result = subprocess.run(
        ["git", "branch", "--show-current"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def indent_block(block: str, indent: str) -> str:
    return "\n".join(indent + line if line else "" for line in block.splitlines())


def update(path: Path, transform: Callable[[str], str]) -> bool:
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    original = path.read_text(encoding="utf-8")
    modified = transform(original)
    if modified == original:
        return False
    path.write_text(modified, encoding="utf-8")
    return True


SEDIMENT_LAYER_FUNCTION = r'''"""
    sediment_layer(deposition, disintegration, time_index;
                   depth=0.0, thickness=1.0,
                   amount_to_cells=Float64)

Reconstruct the sediment stack through `time_index` and return a preserved
interval below the sediment surface.

`deposition` and `disintegration` have dimensions `(facies, x, y, time)`.
`depth` and `thickness` are expressed in sediment-buffer cells.
`amount_to_cells` converts one sediment amount to the same dimensionless units.
The result is `(layer, present)`, where `layer` has dimensions
`(facies, x, y)` and `present` marks cells that contain the requested interval.
"""
function sediment_layer(
    deposition::AbstractArray{T,4},
    disintegration::AbstractArray{T,4},
    time_index::Integer;
    depth::Real = 0.0,
    thickness::Real = 1.0,
    amount_to_cells = Float64,
) where T
    size(deposition) == size(disintegration) ||
        throw(DimensionMismatch(
            "deposition and disintegration must have the same shape",
        ))

    n_facies, nx, ny, n_times = size(deposition)
    1 <= time_index <= n_times ||
        throw(ArgumentError("time_index must be between 1 and $(n_times)"))
    depth >= 0.0 || throw(ArgumentError("depth must be non-negative"))
    thickness > 0.0 || throw(ArgumentError("thickness must be positive"))

    function parcel_mass(data, i, j, k)
        mass = 0.0
        @inbounds for f in 1:n_facies
            amount = Float64(amount_to_cells(data[f, i, j, k]))
            amount >= 0.0 ||
                throw(ArgumentError("sediment amounts must be non-negative"))
            mass += amount
        end
        return mass
    end

    # Determine only the stack depth needed to preserve the requested final
    # interval. In the absence of erosion this is approximately
    # depth + thickness, rather than the complete accumulated succession.
    required_capacity = depth + thickness
    @inbounds for j in 1:ny, i in 1:nx
        height = 0.0
        maximum_height = 0.0

        for k in 1:time_index
            height = max(
                0.0,
                height - parcel_mass(disintegration, i, j, k),
            )
            height += parcel_mass(deposition, i, j, k)
            maximum_height = max(maximum_height, height)
        end

        if height > depth
            sample_bottom = max(0.0, height - depth - thickness)
            required_capacity = max(
                required_capacity,
                maximum_height - sample_bottom,
            )
        end
    end

    # Extra empty rows keep complete-buffer pops away from the edge case in
    # pop_sediment! while leaving the reconstructed result unchanged.
    n_layers = max(3, ceil(Int, required_capacity) + 3)

    layer = zeros(Float64, n_facies, nx, ny)
    present = falses(nx, ny)
    column = zeros(Float64, n_layers, n_facies)
    parcel = zeros(Float64, n_facies)

    @inbounds for j in 1:ny, i in 1:nx
        fill!(column, 0.0)
        available = 0.0

        for k in 1:time_index
            eroded = min(
                parcel_mass(disintegration, i, j, k),
                available,
            )

            if eroded > 0.0
                pop_sediment!(column, eroded)
                available -= eroded
            end

            deposited = 0.0
            for f in 1:n_facies
                amount = Float64(amount_to_cells(deposition[f, i, j, k]))
                amount >= 0.0 ||
                    throw(ArgumentError("sediment amounts must be non-negative"))
                parcel[f] = amount
                deposited += amount
            end

            if deposited > 0.0
                push_sediment!(column, parcel)
                available = min(Float64(n_layers), available + deposited)
            end
        end

        available <= depth && continue

        if depth > 0.0
            pop_sediment!(column, depth)
            available -= depth
        end

        sampled = min(thickness, available)
        sampled <= 0.0 && continue

        layer[:, i, j] .= pop_sediment!(column, sampled)
        present[i, j] = true
    end

    return layer, present
end'''


SEDIMENT_LAYER_TEST = r'''@testset "Sediment layer reconstruction" begin
    using CarboKitten.SedimentStack: sediment_layer

    deposition = zeros(Float64, 3, 1, 1, 4)
    disintegration = zeros(Float64, 3, 1, 1, 4)

    deposition[1, 1, 1, 1] = 1.0
    deposition[2, 1, 1, 2] = 1.0
    deposition[3, 1, 1, 3] = 1.0

    top, present = sediment_layer(deposition, disintegration, 3)
    @test present[1, 1]
    @test top[:, 1, 1] ≈ [0.0, 0.0, 1.0]

    middle, present = sediment_layer(deposition, disintegration, 3; depth=1.0)
    @test present[1, 1]
    @test middle[:, 1, 1] ≈ [0.0, 1.0, 0.0]

    bottom, present = sediment_layer(deposition, disintegration, 3; depth=2.0)
    @test present[1, 1]
    @test bottom[:, 1, 1] ≈ [1.0, 0.0, 0.0]

    combined, present = sediment_layer(
        deposition,
        disintegration,
        3;
        thickness=2.0,
    )
    @test present[1, 1]
    @test combined[:, 1, 1] ≈ [0.0, 1.0, 1.0]

    earlier, present = sediment_layer(deposition, disintegration, 2)
    @test present[1, 1]
    @test earlier[:, 1, 1] ≈ [0.0, 1.0, 0.0]

    absent, present = sediment_layer(deposition, disintegration, 3; depth=3.0)
    @test !present[1, 1]
    @test iszero(sum(absent[:, 1, 1]))

    eroded_deposition = zeros(Float64, 2, 1, 1, 3)
    eroded_disintegration = zeros(Float64, 2, 1, 1, 3)
    eroded_deposition[1, 1, 1, 1] = 1.0
    eroded_deposition[2, 1, 1, 2] = 1.0
    eroded_disintegration[2, 1, 1, 3] = 1.0

    exposed, present = sediment_layer(
        eroded_deposition,
        eroded_disintegration,
        3,
    )
    @test present[1, 1]
    @test exposed[:, 1, 1] ≈ [1.0, 0.0]

    # A finite layer thickness combines a very small final event with the
    # underlying sediment instead of letting that event dominate the map.
    thin_deposition = zeros(Float64, 2, 1, 1, 2)
    thin_disintegration = zeros(Float64, 2, 1, 1, 2)
    thin_deposition[1, 1, 1, 1] = 1.0
    thin_deposition[2, 1, 1, 2] = 0.01

    smoothed, present = sediment_layer(
        thin_deposition,
        thin_disintegration,
        2;
        thickness=1.0,
    )
    @test present[1, 1]
    @test smoothed[:, 1, 1] ≈ [0.99, 0.01]
end'''


FACIES_FRACTION_FUNCTION = r'''function _facies_fraction(d::AbstractArray, facies::Integer)
    total = dropdims(sum(d; dims = 1); dims = 1)
    selected = d[facies, :, :]
    return ifelse.(
        iszero.(total),
        missing,
        Float64.(ustrip.(selected ./ total)),
    )
end'''


PRESERVED_VALUES_FUNCTION = r'''function preserved_values()
    layer, present = sediment_layer(
        data.deposition,
        data.disintegration,
        t_idx;
        depth = depth_cells,
        thickness = layer_thickness_cells,
        amount_to_cells = amount_to_cells,
    )

    m = if color_by == :facies
        Matrix{Union{Missing,Int}}(_colormax(layer))
    else
        _facies_fraction(layer, facies_idx)
    end

    m[.!present] .= missing
    if mask_emerged
        m[wd .> 0u"m"] .= missing
    end

    return m
end'''


def patch_sediment_module(text: str, label: str) -> str:
    if "function sediment_layer(" in text:
        return text

    export_pattern = re.compile(
        r"(?m)^(?P<indent>[ \t]*)export push_sediment!, pop_sediment!, peek_sediment\s*$"
    )
    matches = list(export_pattern.finditer(text))
    if len(matches) != 1:
        fail(f"{label}: expected one SedimentStack export line, found {len(matches)}")

    match = matches[0]
    replacement = (
        match.group("indent")
        + "export push_sediment!, pop_sediment!, peek_sediment, sediment_layer"
    )
    text = text[: match.start()] + replacement + text[match.end() :]

    marker = re.search(
        r"(?m)^(?P<indent>[ \t]*)function pop_sediment!\("
        r"cols::AbstractArray\{F,\s*4\}",
        text,
    )
    if marker is None:
        fail(f"{label}: could not find the 4-D pop_sediment! method")

    block = indent_block(SEDIMENT_LAYER_FUNCTION, marker.group("indent")) + "\n\n"
    return text[: marker.start()] + block + text[marker.start() :]


def patch_sediment_test_file(text: str, label: str) -> str:
    if '@testset "Sediment layer reconstruction"' in text:
        return text

    end_marker = re.search(r"(?m)^# ~/~ end\s*$", text)
    if end_marker is None:
        insertion = len(text.rstrip())
        suffix = "\n"
    else:
        insertion = end_marker.start()
        suffix = ""

    prefix = text[:insertion].rstrip() + "\n\n"
    return prefix + SEDIMENT_LAYER_TEST + "\n" + suffix + text[insertion:]


def patch_sediment_doc(text: str) -> str:
    text = patch_sediment_module(text, str(SEDIMENT_DOC))

    if '@testset "Sediment layer reconstruction"' not in text:
        heading = text.find("### Pushing sediment")
        if heading < 0:
            fail("sediment-buffer documentation: missing '### Pushing sediment'")

        test_start = text.rfind('@testset "SedimentArray"', 0, heading)
        if test_start < 0:
            fail("sediment-buffer documentation: missing SedimentArray test")

        line_start = text.rfind("\n", 0, test_start) + 1
        indent = text[line_start:test_start]

        tail = text[test_start:heading]
        test_line = re.search(
            r"(?m)^[ \t]*@test all\(sum\(a; dims=1\) \.≈ 1\.0\)\s*$",
            tail,
        )
        if test_line is None:
            fail("sediment-buffer documentation: missing SedimentArray assertion")

        after_assertion = test_start + test_line.end()
        end_match = re.search(
            r"(?m)^" + re.escape(indent) + r"end\s*$",
            text[after_assertion:heading],
        )
        if end_match is None:
            fail("sediment-buffer documentation: missing SedimentArray test end")

        insertion = after_assertion + end_match.end()
        block = "\n" + indent_block(SEDIMENT_LAYER_TEST, indent)
        text = text[:insertion] + block + text[insertion:]

    if "Map views reconstruct stratigraphy" not in text:
        heading = "## Component"
        position = text.find(heading)
        if position < 0:
            fail("sediment-buffer documentation: missing Component heading")

        explanation = '''### Extracting a stratigraphic layer\n\nMap views reconstruct stratigraphy by replaying deposition and disintegration\nthrough the same sediment-stack operations used by the model. After the stack\nhas been built to the requested time, `sediment_layer` removes the selected\noverburden and returns a finite interval. The finite thickness prevents a very\nsmall last sedimentation event from producing unstable facies fractions.\n\nThe helper is unit-free, like the rest of `SedimentStack`. A caller working in\nphysical lengths converts sediment amounts, depth, and layer thickness with the\nmodel's depositional resolution.\n\n'''
        text = text[:position] + explanation + text[position:]

    return text


def replace_indented_function(
    text: str,
    name_pattern: str,
    replacement: str,
    following_pattern: str,
    label: str,
) -> str:
    pattern = re.compile(
        rf"(?ms)^(?P<indent>[ \t]*)function {name_pattern}.*?"
        rf"^(?P=indent)end\s*\n(?=(?P<following>[ \t]*{following_pattern}))"
    )
    match = pattern.search(text)
    if match is None:
        fail(f"{label}: could not find function {name_pattern}")
    block = indent_block(replacement, match.group("indent")) + "\n"
    return text[: match.start()] + block + text[match.end() :]


def patch_map_source(text: str, label: str) -> str:
    if "using CarboKitten.SedimentStack: sediment_layer" not in text:
        pattern = re.compile(
            r"(?m)^(?P<indent>[ \t]*)using CarboKitten\.Output\.Abstract: "
            r"stratigraphic_column, water_depth\s*$"
        )
        match = pattern.search(text)
        if match is None:
            fail(f"{label}: could not find the Output.Abstract import")
        indent = match.group("indent")
        replacement = (
            indent
            + "using CarboKitten.Output.Abstract: water_depth\n"
            + indent
            + "using CarboKitten.SedimentStack: sediment_layer"
        )
        text = text[: match.start()] + replacement + text[match.end() :]

    if "return ifelse.(" not in text:
        text = replace_indented_function(
            text,
            r"_facies_fraction\(d::AbstractArray, facies::Integer\)",
            FACIES_FRACTION_FUNCTION,
            r"# Resolve a stratigraphic position",
            label,
        )

    # Add the new keywords to the displayed signature in the docstring.
    if re.search(r"(?m)^[ \t]*depth\s*=\s*0\.0u\"m\",\s*$", text) is None:
        doc_time = re.search(
            r"(?m)^(?P<indent>[ \t]*)time(?P<spacing>[ \t]*)=[ \t]*end,\s*$",
            text,
        )
        if doc_time is None:
            fail(f"{label}: could not find docstring time keyword")
        indent = doc_time.group("indent")
        spacing = doc_time.group("spacing")
        width = max(len("time") + len(spacing), len("depositional_resolution"))

        def aligned(name: str, value: str) -> str:
            return indent + name.ljust(width) + " = " + value + ","

        addition = "\n".join(
            [
                aligned("depth", '0.0u"m"'),
                aligned("layer_thickness", '1.0u"m"'),
                aligned("depositional_resolution", '0.5u"m"'),
            ]
        )
        text = text[: doc_time.end()] + "\n" + addition + text[doc_time.end() :]

    if "depth = 0.0u\"m\"," not in text:
        signature = re.search(
            r"(?m)^(?P<indent>[ \t]*)time::Union\{Integer,Quantity\}"
            r"[ \t]*=[ \t]*length\(header\.axes\.t\[1:data\.write_interval:end\]\),\s*$",
            text,
        )
        if signature is None:
            fail(f"{label}: could not find map_view! time argument")
        indent = signature.group("indent")
        addition = (
            f'\n{indent}depth = 0.0u"m",'
            f'\n{indent}layer_thickness = 1.0u"m",'
            f'\n{indent}depositional_resolution = 0.5u"m",'
        )
        text = text[: signature.end()] + "".join(addition) + text[signature.end() :]

    if "depth_cells = Float64(ustrip(depth_m / resolution_m))" not in text:
        prec = re.search(
            r"(?m)^(?P<indent>[ \t]*)prec = 1e-8u\"m\".*$",
            text,
        )
        if prec is None:
            fail(f"{label}: could not find the old precision line")
        indent = prec.group("indent")
        validation = r'''depth_m = uconvert(u"m", depth)
layer_thickness_m = uconvert(u"m", layer_thickness)
resolution_m = uconvert(u"m", depositional_resolution)

depth_m >= 0.0u"m" || error("`depth` must be non-negative")
layer_thickness_m > 0.0u"m" ||
    error("`layer_thickness` must be positive")
resolution_m > 0.0u"m" ||
    error("`depositional_resolution` must be positive")

depth_cells = Float64(ustrip(depth_m / resolution_m))
layer_thickness_cells = Float64(ustrip(layer_thickness_m / resolution_m))
amount_to_cells = amount -> Float64(ustrip(amount / resolution_m))'''
        block = indent_block(validation, indent)
        text = text[: prec.start()] + block + text[prec.end() :]

    if "layer, present = sediment_layer(" not in text:
        text = replace_indented_function(
            text,
            r"preserved_values\(\)",
            PRESERVED_VALUES_FUNCTION,
            r"# Merge defaults",
            label,
        )

    if "depth_value = ustrip(depth_m)" not in text:
        title = re.search(
            r"(?m)^(?P<indent>[ \t]*)ax\.title[ \t]*=[ \t]*"
            r'"t = \$\(round\(t_myr; digits = 3\)\) Myr"\s*$',
            text,
        )
        if title is None:
            fail(f"{label}: could not find the map title")
        indent = title.group("indent")
        title_block = r'''if show == :model
    ax.title = "t = $(round(t_myr; digits = 3)) Myr"
else
    depth_value = ustrip(depth_m)
    ax.title = "t = $(round(t_myr; digits = 3)) Myr, " *
               "depth = $(round(depth_value; digits = 3)) m"
end'''
        block = indent_block(title_block, indent)
        text = text[: title.start()] + block + text[title.end() :]

    return text


def strip_entangled_markers(source: str) -> str:
    lines = source.splitlines()
    while lines and lines[0].startswith("# ~/~ begin"):
        lines.pop(0)
    while lines and lines[-1].startswith("# ~/~ end"):
        lines.pop()
    return "\n".join(lines).strip("\n")


def sync_map_implementation(doc: str, source: str) -> str:
    heading = doc.find("### Implementation")
    if heading < 0:
        fail("map-view documentation: missing Implementation heading")

    module_match = re.search(
        r"(?m)^(?P<indent>[ \t]*)module MapView\s*$",
        doc[heading:],
    )
    if module_match is None:
        fail("map-view documentation: missing MapView module block")

    start = heading + module_match.start()
    indent = module_match.group("indent")
    end_pattern = re.compile(
        r"(?m)^" + re.escape(indent) + r"end\s*#\s*module MapView\s*$"
    )
    end_match = end_pattern.search(doc, start)
    if end_match is None:
        fail("map-view documentation: missing end of MapView module block")

    implementation = indent_block(strip_entangled_markers(source), indent)
    return doc[:start] + implementation + doc[end_match.end() :]


def patch_map_doc_prose(text: str) -> str:
    if "replays all saved deposition and disintegration" not in text:
        anchor = (
            "Maps can be coloured either categorically, by dominant facies, or "
            "continuously, by the proportion of one selected facies relative to "
            "the total sediment in each cell."
        )
        if anchor not in text:
            fail("map-view documentation: missing introduction paragraph")
        addition = '''\n\nFor preserved maps, the routine replays all saved deposition and disintegration\nthrough the existing sediment-buffer operations up to the selected time. It\nthen samples a finite interval at the requested depth below the sediment\nsurface. This separates stratigraphic reconstruction from plotting and prevents\na very small last sedimentation event from dominating facies fractions.\n\nFor exact timestep-by-timestep reconstruction, use an output\n`write_interval = 1`. Larger intervals remain supported because CarboKitten\naggregates deposition and disintegration within each saved interval.'''
        text = text.replace(anchor, anchor + addition, 1)

    if "- `depth` — depth below the sediment surface" not in text:
        pattern = re.compile(
            r"(?ms)^(?P<indent>[ \t]*)- `time` —.*?"
            r"the nearest available frame is used\.[ \t]*$"
        )
        match = pattern.search(text)
        if match is None:
            fail("map-view documentation: missing time-keyword description")
        indent = match.group("indent")
        bullets = '''- `depth` — depth below the sediment surface at the selected time. The
  default is `0u"m"`, i.e. the uppermost preserved sediment.
- `layer_thickness` — thickness of the sampled interval. The default is
  `1u"m"`, which reduces noise from very small sedimentation events.
- `depositional_resolution` — vertical resolution used for the reconstructed
  sediment buffer. The default is `0.5u"m"`; use the value from the model input
  when a different resolution was used.'''
        block = "\n" + indent_block(bullets, indent)
        text = text[: match.end()] + block + text[match.end() :]

    return text


def run_tangle() -> bool:
    commands: list[list[str]] = []
    if shutil.which("uv") is not None:
        commands.append(["uv", "run", "entangled", "tangle"])
    if shutil.which("entangled") is not None:
        commands.append(["entangled", "tangle"])

    for command in commands:
        print("Tangling generated files:", " ".join(command))
        result = subprocess.run(command, check=False)
        if result.returncode == 0:
            return True

    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force-branch",
        action="store_true",
        help="apply even if the current branch name differs",
    )
    parser.add_argument(
        "--no-tangle",
        action="store_true",
        help="edit only the Entangled Markdown sources",
    )
    args = parser.parse_args()

    branch = current_branch()
    if branch != BRANCH and not args.force_branch:
        fail(
            f"Current branch is {branch!r}. Switch to {BRANCH!r}, "
            "or rerun with --force-branch."
        )

    changed: list[str] = []

    if update(SEDIMENT_DOC, patch_sediment_doc):
        changed.append(str(SEDIMENT_DOC))

    def transform_map_doc(text: str) -> str:
        text = patch_map_doc_prose(text)
        return patch_map_source(text, str(MAP_DOC))

    if update(MAP_DOC, transform_map_doc):
        changed.append(str(MAP_DOC))

    required_docs = {
        SEDIMENT_DOC: "function sediment_layer(",
        MAP_DOC: "layer, present = sediment_layer(",
    }
    for path, marker in required_docs.items():
        if marker not in path.read_text(encoding="utf-8"):
            fail(f"Verification failed: {marker!r} missing from {path}")

    tangled = False
    if not args.no_tangle:
        tangled = run_tangle()
        if not tangled:
            print(
                "\nEntangled was not available or tangling failed. "
                "The Markdown sources are updated; run one of:\n"
                "  uv run entangled tangle\n"
                "  entangled tangle"
            )

    if tangled:
        required_generated = {
            SEDIMENT_SRC: "function sediment_layer(",
            SEDIMENT_TEST: '@testset "Sediment layer reconstruction"',
            MAP_SRC: "layer, present = sediment_layer(",
        }
        for path, marker in required_generated.items():
            if not path.is_file() or marker not in path.read_text(encoding="utf-8"):
                fail(f"Tangling verification failed: {marker!r} missing from {path}")

    print("\nUpdated files:" if changed else "\nPatch already applied; no Markdown files changed.")
    for path in changed:
        print(f"  {path}")

    print("\nRun next:")
    if args.no_tangle:
        print("  uv run entangled tangle")
    print("  git diff --check")
    print("  julia --project=. -e 'using Pkg; Pkg.test()'")
    print(
        "  git diff -- docs/src/components/sediment_buffer.md "
        "src/SedimentStack.jl test/SedimentStackSpec.jl "
        "docs/src/visualization/map-view.md ext/MapView.jl"
    )


if __name__ == "__main__":
    main()
