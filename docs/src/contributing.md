# Contributing

- Try to explain your contribution with simple language.
- References are always welcome.
- Follow the coding standards in the source.

## Reporting issues

If you are experiencing issues or have discovered a bug, please [report it on GitHub](https://github.com/henry2004y/Vlasiator.jl/issues).
To make the resolution process easier, please include the version of Julia and Vlasiator.jl in your writeup.
These can be found with two commands:

```julia
julia> versioninfo()
julia> using Pkg; Pkg.status()
```

## Feature requests

If you have suggestions of improvement or algorithms that you would like to see implemented in Vlasiator.jl, please open an [issue](https://github.com/henry2004y/Vlasiator.jl/issues) on GitHub. Suggestions as well as feature requests are very welcome.

## Code contribution

If you have code that you would like to contribute to Vlasiator.jl, that is awesome! Please open an [issue](https://github.com/henry2004y/Vlasiator.jl/issues)
before you create the pull request on GitHub so that we make sure your idea is aligned with our goals for the project.

After your idea is discussed and revised by maintainers, please get the development version of the project by typing the following in the package manager:

```julia
] activate @dev
```

This will create a fresh environment called `@dev` where you can play with the project components without compromising your normal
user environment.

```julia
] dev Vlasiator
```

This will clone all the project components in your `~/.julia` folder so that you can modify it and submit a pull request on GitHub later. Don't hesitate to ask questions. We are looking forward to your contributions.