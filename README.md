# HighOrderMeshes

[![Build Status](https://github.com/popersson/HighOrderMeshes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/popersson/HighOrderMeshes.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Installation and Versions

Since this package is in active development, the `main` branch may contain breaking changes.

### For old stable version (2024-2025 Legacy Code)

If you are working on a project that uses the original interface (from 2024-2025), you should use the stable **v0.1.0** tag. To switch to this version, run the following commands in your terminal inside the project folder:

```bash
# Switch to the stable v0.1.0 version
git checkout v0.1.0

# Update the Julia environment
julia --project -e 'using Pkg; Pkg.instantiate()'

```

> **Note:** When you run the checkout command, Git will mention a "detached HEAD" state. This is normal and simply means you are locked to that specific point in time.

### Current Development

To stay on the latest version with all the newest features, ensure you are on the `main` branch:

```bash
git checkout main

```
