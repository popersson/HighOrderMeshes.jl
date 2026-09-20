using Documenter, HighOrderMeshes

makedocs(
    sitename = "HighOrderMeshes.jl",
    modules  = [HighOrderMeshes],
    format   = Documenter.HTML(
        edit_link = "main",
        canonical = "https://popersson.github.io/HighOrderMeshes.jl/stable",
    ),
    pages    = [
        "Home"          => "index.md",
        "Quick start"   => "quickstart.md",
        "Mesh format"   => "meshformat.md",
        "Basic meshes"  => "basicmeshes.md",
        "Solver node sets" => "solvernodes.md",
        "Gmsh import"   => "gmsh.md",
        "File I/O"      => "io.md",
        "API reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/popersson/HighOrderMeshes.jl.git",
    devbranch = "main",
    push_preview = true,
)
