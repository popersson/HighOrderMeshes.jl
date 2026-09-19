# HighOrderMeshes.jl

![CI](https://github.com/popersson/HighOrderMeshes.jl/actions/workflows/CI.yml/badge.svg)
[![codecov](https://codecov.io/gh/popersson/HighOrderMeshes.jl/graph/badge.svg?token=FLXZ69IRUK)](https://codecov.io/gh/popersson/HighOrderMeshes.jl)
[![docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://popersson.github.io/HighOrderMeshes.jl/stable)
[![docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://popersson.github.io/HighOrderMeshes.jl/dev)
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

### Earlier versions

Version 0.3.0 is a breaking redesign of the core API; see `CHANGELOG.md`.
The last version before the redesign is `v0.2.0`, and the original 2024-2025
interface is `v0.1.0`. Install either by specifying the tag:

```julia
Pkg.add(url="https://github.com/popersson/HighOrderMeshes.jl", rev="v0.2.0")
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
