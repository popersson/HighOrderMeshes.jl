# HighOrderMeshes.jl

![CI](https://github.com/popersson/HighOrderMeshes.jl/actions/workflows/CI.yml/badge.svg)
Tools for high-order unstructured meshes and finite element methods.
[![codecov](https://codecov.io/gh/popersson/HighOrderMeshes.jl/graph/badge.svg?token=FLXZ69IRUK)](https://codecov.io/gh/popersson/HighOrderMeshes.jl)

## Installation

This package is not yet in the General Registry. You can install the latest version directly from GitHub:

```julia
import Pkg
Pkg.add(url="https://github.com/popersson/HighOrderMeshes.jl")
```

> **Warning:** Since this package is in active development, the `main` branch may contain breaking changes.

### Legacy Version (2024-2025 Code)

If you need the original interface (v0.1.0), install it by specifying the tag:

```julia
Pkg.add(url="https://github.com/popersson/HighOrderMeshes.jl", rev="v0.1.0")
```

## Quick Start

### Visualization (Makie.jl)

We recommend **Makie.jl** for visualization.

```julia
using HighOrderMeshes
using GLMakie # Or CairoMakie for non-interactive plots

msh = ex1mesh()
plot(msh)           # Plot high-order mesh
```

```julia
u = ex1solution(msh)
plot(msh, u)        # Plot sample solution
```

### Visualization (Plots.jl)

If you prefer Plots.jl, you must install and load `TriplotRecipes.jl` for the plotting extension to activate.

```julia
using HighOrderMeshes
using Plots, TriplotRecipes

plot(ex1mesh())

```
