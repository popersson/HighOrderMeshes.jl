using Documenter, HighOrderMeshes

makedocs(
    sitename = "HighOrderMeshes.jl",
    modules  = [HighOrderMeshes],
    pages    = [
        "Home"          => "index.md",
        "Quick start"   => "quickstart.md",
        "Mesh format"   => "meshformat.md",
        "Basic meshes"  => "basicmeshes.md",
        "Gmsh import"   => "gmsh.md",
        "File I/O"      => "io.md",
        "API reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/popersson/HighOrderMeshes.jl.git",
    push_preview = true,
)
