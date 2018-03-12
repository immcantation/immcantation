We welcome contributions to all components of the Immcantation framework through
pull request to the relevant Bitbucket repository:

+ [Docker container, pipelines and portal documentation](https://bitbucket.org/kleinstein/immcantation)
+ [pRESTO](https://bitbucket.org/kleinstein/presto)
+ [Change-O](https://bitbucket.org/kleinstein/changeo)
+ [Alakazam](https://bitbucket.org/kleinstein/alakazam)
+ [SHazaM](https://bitbucket.org/kleinstein/shazam)
+ [TIgGER](https://bitbucket.org/kleinstein/tigger)

All packages are under the free and open-source license
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0).

---

# Documentation Guidelines

---

### Python Package Documentation

+ Use [reStructuredText](http://www.sphinx-doc.org/en/stable/rest.html) to write documentation.
+ Use [Sphinx](http://www.sphinx-doc.org/en/stable) to build documentation.
+ Document functions with [Google style](https://www.chromium.org/chromium-os/python-style-guidelines)
  docstrings, except use 2 space indents and spell out 'Arguments' instead of using 'Args'.
+ Documentation for commandline arguments are generating using [argparse](https://docs.python.org/3/library/argparse.html),
  which requires all commandline arguments to be defined in a function that returns the `argparse.ArgumentParser` object.
+ API documentation is generated using [autodoc](http://www.sphinx-doc.org/en/stable/ext/autodoc.html) with
  [napoleon](http://www.sphinx-doc.org/en/stable/ext/napoleon.html).
+ The package [autoprogram](http://pythonhosted.org/sphinxcontrib-autoprogram) is used to automatically generate commandline
  documentation. This requires each Python script to define all commandline arguments in a function that returns an
  [argparse.ArgumentParser](https://docs.python.org/3/library/argparse.html) object.
+ The README.rst, INSTALL.rst and NEWS.rst are automatically include in the docs via the Sphinx
  `include` directive (eg, `.. include:: ../README.rst`), so changes to this content should be make in these top level files.

### R Package Documentation

+ Write vignettes in [R Markdown](http://rmarkdown.rstudio.com).
+ Build vignettes using [knitr](http://yihui.name/knitr).
+ Document functions, classes, methods and data with [Roxygen](http://r-pkgs.had.co.nz/man.html) tags which will
  automatically generate man pages from the function documentation.
+ Use [markr](https://bitbucket.org/javh/markr) to convert the R man pages, vignettes, DESCRIPTION, README and
  CITATION files to Markdown to build the [MkDocs](http://www.mkdocs.org) documentation that is hosted on ReadTheDocs.

### Best practices for rebuilding R docs with markr for an existing package:

1. Make sure the date, version and changelog in DESCRIPTION and NEWS are current!
2. Delete the contents of `docs/topics` and `docs/vignettes`
3. Run the script in `inst/markr/build.R`.
4. Update `mkdocs.yml` with any topics/vignettes that have been added/removed.
5. Verify the docs build correctly by running `mkdocs serve` from the package directory and going to the url specified.
6. Commit and push. Making sure to add/remove files from the repo that are new/lost.
7. Verify the ReadTheDocs site builds correctly.

---

# Unit Test Guidelines

---

### Python Unit Tests

+ Use the [unittest](https://docs.python.org/3/library/unittest.html) framework.
+ Create separate files for each module inside the `tests` folder.
+ If data is required for a unit test, place it in `tests/data` and access the
  directory via the recipe:
  ```
  test_path = os.path.dirname(os.path.realpath(__file__))
  data_path = os.path.join(test_path, 'data')
  ```

### R Unit Tests

+ Use the [testthat](https://github.com/hadley/testthat) framework.
+ Place tests in `tests/testthat` within the R file corresponding
  to the file in `/R` that contains the function being tested.
+ If data is required for a unit test, place it in `tests/data-tests`
  and access it via the recipe:
  ```
  data <- file.path("..", "data-tests", "data.rda")
  ```

---

# Coding Style Guidelines

---

We tend to use the same style for both R and Python. This style differs
from the canonical styles for both Python and Bioconductor, but is largely
similar to both. It is essentially a compromise between these two
standards and very close to what you typically see for Java and C++.

### Indentation

+ Use 4 spaces for indenting. No tabs.
+ No lines longer than 90 characters.  Set the line marker in RStudio
  and PyCharm accordingly. Consider 80 characters a soft limit. Try
  to stick to less than 80 within code, but let it stretch to 90 in
  roxygen blocks, docstrings and argparse help if you have to.
  Readability should be the priority.

### Variable names

+ Use lower case with underscore separation: `variable_name`.
+ For constants and variables in shell scripts use upper case with
  underscore separation: `CONSTANT_NAME`, `SHELL_VARIABLE`.
+ It's often convenient to specify the variable type in the name if you
  are transforming the data (Hungarian notation).  For example,
  converting the data.frame `data_df` to the list `data_list`, without
  much else happening. But don't overdo it. Simple is better if
  the type and purpose of the variable is obvious.
+ Concise, but meaningful, variable names are preferred, except for
  variables in list comprehensions, for loop indices, lambda functions,
  etc, where variables like `x` or `i` would be cleaner.
+ Single word names for function arguments are preferable, if it's
  possible to do so without being confusing or colliding with reserved
  words and existing function/builtin names.

### Function and method names

+ Use lower camelCase: `functionName`.
+ Use verbNoun for function names (eg, `getGene`) if it makes sense to do so.
+ Do not use periods in function names within R. In the S3 class system,
  `plot(x)`, where `x` is class `Foo`, will dispatch to `plot.Foo(x)`.

### Classes, R data, Python module names

+ Use upper CamelCase: `ClassName`, `DataObject`, `PythonModule`.

### Spaces

+ Always use space after a comma, semicolon, pipe or other delimiter.
  For example: `c(x, y, z)`, `[, 5]`, `cat foo.txt | grep bar`.
+ No space around `=` when using named arguments to functions:
  `functionName(x=1, y=2)`
+ Space around most binary operators: `x == y`, `x + y`. However, it is
  often more readable to have high priority operators without a space and
  low priority with: `2*x + y^4`.
+ Don't put spaces within the inner edge of parenthesis:
  `functionName(x, y)`, `(z - (x + y))`
+ Don't put spaces between a function and its parenthesis: `functionName(x, y)`
+ Put a space around expressions within parenthesis:
  `if (x + y == 4) { print(TRUE) }`.

### Braces in R

+ The opening brace should be on the same line as the function declaration or
  for/if/while statement, and the closing brace should be on a line of its own
  indented to match the position of the statement.
+ Braces should always be used for function declarations and for/if/while
  statements even if they are only one line: `for (x in 1:10) { print(x) }`

### Assignment in R

+ Use `<-`, not `=`, for variable assignment.
+ Use `=` for assigning defaults in function declarations.

### Comments

+ Indent at the same level as surrounding code.
+ Fully document each function with its purpose, arguments and return
  values in either the docstring or roxygen block.
+ Comments within code should act sort of like a table of contents. They
  should tell you what the purpose of the code is, but not regurgitate
  the implementation. They are also good for acting like section headers.
  Basically, you want to be able to skim the comments and understand what
  the code does and find a specific section you are looking for without
  having to read the code itself.
+ Comments should also be used to explain tricky or confusing lines of
  code, but in general it is better to avoid writing this type of code
  in the first place.
+ Use the RStudio section header syntax (`#### Section name ####`) to
  separate R code into logical sections.

### Everything else

For everything else see either the
[Bioconductor](http://bioconductor.org/developers/how-to/coding-style)
or [Python](https://www.python.org/dev/peps/pep-0008) guides.

---

# Version Control Guidelines

---

### Hidden Files ###

As a general rule, do not commit hidden (`.*`) files such as project build files
(eg, `*.Rproj`, `.idea`) or other local environment settings (eg, `.hgignore`),
with the noted exception of `.Rbuildignore`. We use the following
standard `.hgignore` files locally:

Standard `.hgignore` file for the Python packages:

```
syntax: glob
.*
*~
*.pyc
*.pyo
*.egg-info
run
run_*
dist
build
docs/_build
```

Standard `.hgignore` file for the R packages:

```
syntax: glob
.*
*.Rproj
man/*.Rd
inst/doc
vignettes/*.R
vignettes/*.html
vignettes/*.pdf
src/*.o
src/*.so
src/*.dll
```

### Branches ###

* The `default` branch acts as the main development branch. Try to avoid putting
  non-working code in the `default` branch.
* Used named branches for adding large features or for work that will distrupt the stability of
  the `default` branch for more than a day or two.
* If necessary, open a maintenance branch to update an old, backwards incompatible
  version major version. For example, `v0.3_maintenance` to update code in v0.3 when default
  has moved on to v0.4.

### Version Numbers

+ Denote the version number and date in the `package/Version.py` file for Python packages
  and the `DESCRIPTION` file for R packages.
+ Use the `x.y.z` style of version numbers, usually starting with 0.1.0, and following
  the convention `<major>.<minor>.<patch>`. An increase in `<major>` usually means
  large changes and/or breaking of backwards compatibility, `<patch>` usually means
  bug fixes and small changes, with `<minor>` falling somewhere in-between.
+ Label development code by updating the version in the package to `x.y.z.999`, where `x.y.z`
  is the last official release. Do this with the first commit after an official release.
  Add an appropriate development build entry to the NEWS file.
+ Make sure you keep the date in the package current! Because we use a simple versioning
  scheme that doesn't distinguish development builds, the date is your only indication
  of which development build you are using.

### Annotating Versions

* Use tags in mercurial to denote official versions in the form `Version 0.1.1 - Description`.
* Do not tag development builds. Just leave `x.y.z.999` in the source code.
* Before tagging an official release:
    1. Check that all easily forgotten files are up to date, including README, INSTALL, NEWS,
       DESCRIPTION, and so on.
    2. Make sure the package documentation is all current, including building and validating docs for
       [ReadTheDocs](http://readthedocs.io).
    3. Do not tag the release until **after** the package has been accepted by CRAN or PyPI.
* After an official version has been tagged, head over to the ReadTheDocs project,
  go to Admin | Settings | Versions, and activate the version.