# HighOrderMeshes.jl

![CI](https://github.com/popersson/HighOrderMeshes.jl/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/popersson/HighOrderMeshes.jl/graph/badge.svg?token=FLXZ69IRUK)](https://codecov.io/gh/popersson/HighOrderMeshes.jl)
[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://popersson.github.io/HighOrderMeshes.jl)
Tools for high-order unstructured meshes and finite element methods.

See the [documentation](https://popersson.github.io/HighOrderMeshes.jl) for
a full quick start, the mesh format, Gmsh import, and the API reference.

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

Visualization is provided through a [Makie.jl](https://docs.makie.org)
package extension; load any Makie backend to enable `plot`.

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
