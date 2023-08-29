using Documenter, CarboKitten

makedocs(
    sitename="CarboKitten",
    pages = [
        "Bosscher and Schlager 1992" => "bosscher-1992.md",
        "CarboCAT" => [
            "summary" => "carbocat.md",
            "cellular automaton" => "carbocat-ca.md",
            "sediment transport" => "carbocat-transport.md"
        ],
        "Algorithms" => [
            "stencils" => "stencils.md",
            "utility" => "utility.md"
        ]
    ])
