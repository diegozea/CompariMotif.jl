```@eval
import CompariMotif
using Markdown
Markdown.parse(read(joinpath(pkgdir(CompariMotif), "README.md"), String))
```

## Public API

```@autodocs
Modules = [CompariMotif]
Private = false
```

## Internal API

```@autodocs
Modules = [CompariMotif]
Public = false
Private = true
```
