using TimerOutputs

dto = TimerOutput()
reset_timer!(dto)

using Ferrite
using FerriteGridUtil
using Documenter

DocMeta.setdocmeta!(FerriteGridUtil, :DocTestSetup, :(using Ferrite, FerriteGridUtil); recursive=true)

# Generate tutorials
include("generate.jl")

makedocs(;
    modules=[FerriteGridUtil],
    authors="David Rollin <d.rollin@tu-braunschweig.de>, Phil Kreutzheide <p.kreutzheide@tu-braunschweig.de>",
    sitename="FerriteGridUtil.jl",
    format=Documenter.HTML(;
        canonical="https://DRollin.github.io/FerriteGridUtil.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "Tutorials overview" => "tutorials/index.md",
            "tutorials/dummy.md",
        ],
        "Reference" => [
            "Reference overview" => "reference/index.md",
            "reference/property.md",
            "reference/manipulation.md",
            "reference/saveload.md",
            "reference/conversion.md",
        ],
        "devdocs/index.md",
    ],
)

deploydocs(;
    repo="github.com/DRollin/FerriteGridUtil.jl.git",
    devbranch="main",
    push_preview=true,
)