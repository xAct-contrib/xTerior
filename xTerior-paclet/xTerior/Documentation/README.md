# xTerior documentation

The editable Wolfram Documentation Tools source notebooks are in:

- `English/ReferencePages/Symbols/`: one page for each of the 39 symbols with an explicit usage message in the kernel source.
- `English/Guides/xTerior.nb`: the main guide, grouped by topic.
- `English/Tutorials/ExteriorCalculus.nb`: **Exterior Calculus with xTerior**, adapted from the supplied `xTeriorDoc.nb` manual.

The pages include usage syntax, implementation notes, See Also links, package initialization, self-contained example setup, evaluated examples, and categorization metadata. The documentation link base is `xTerior`.

## Editing

Open a source notebook in Mathematica and use **Palettes > Documentation Tools**. Edit these source notebooks, not the generated notebooks under `build/`. Run each page's example initialization and setup once in a fresh kernel with xAct installed. xTerior's `Diff` is `xAct\`xTerior\`Diff`, distinct from the recent System symbol of the same name.

## Tutorial

The tutorial preserves the eight main sections of `xTeriorDoc.nb` (version 0.9.1, last updated 19 April 2021), by Alfonso García-Parrado Gómez-Lobo and Leo C. Stein, © 2013–2021 under the GNU General Public License:

1. An example session.
2. Definition and basic properties of differential forms.
3. Cartan structure equations and Chern–Simons forms.
4. The Hodge dual and codifferential.
5. Exterior covariant differentiation.
6. Lie derivatives and interior contraction.
7. Variational derivatives, with electromagnetism, Yang–Mills, and Einstein–Cartan examples.
8. Coordinate calculations, including Schwarzschild geometry and the Hodge Laplacian.

Start a fresh kernel, load xTerior, and evaluate the tutorial's input cells in order. Sections share definitions and use `%` and `%%` to refer to earlier results. The initial example session removes the geometric objects it creates before the subsequent sections introduce their shared setup. The tutorial changes session settings such as the derivation order, dimensional simplification, and canonicalization options; use a dedicated kernel when following the complete session.

The tutorial is linked from the main guide and from the reference pages for the functions it discusses. Its Documentation Center address is `paclet:xTerior/tutorial/ExteriorCalculus`. The original attached manual remains the source attribution; the notebook under `English/Tutorials/` is the maintained paclet documentation source.

## Building and previewing

Use a Mathematica installation with PacletTools, Documentation Tools, and a working notebook front end. xAct is needed to evaluate the examples. In a fresh kernel, set `pacletRoot` to the directory containing `PacletInfo.wl`, then evaluate:

```wolfram
Needs["PacletTools`"];
result = PacletBuild[pacletRoot];
```

`PacletBuild` compiles the documentation and creates `build/xTerior-0.9.1.paclet`. For a documentation-only build, use `PacletDocumentationBuild[pacletRoot]`. Compilation is necessary: simply archiving the editable source notebooks does not produce finished Documentation Center pages.

Preview the built pages without replacing an installed paclet:

```wolfram
PacletDirectoryLoad[FileNameJoin[{pacletRoot, "build", "xTerior"}]];
SystemOpen["paclet:xTerior/guide/xTerior"];
SystemOpen["paclet:xTerior/ref/DefDiffForm"];
SystemOpen["paclet:xTerior/tutorial/ExteriorCalculus"];
```

To install the built archive, evaluate:

```wolfram
PacletInstall[result["PacletArchive"], ForceVersionInstall -> True]
```

`Get[FileNameJoin[{pacletRoot, "makexTerior.wl"}]]` performs the build and installation together. A failed build leaves the installed paclet untouched. The package version remains 0.9.1; increment it before distributing a new release if appropriate.

## Validation and source limitations

Run the notebook coverage and link checks with:

```wolfram
TestReport[FileNameJoin[{pacletRoot, "Tests", "Documentation.wlt"}]]
```

The reference-page examples were originally checked with Wolfram Language 15.0.1, xTensor 1.3.0, xCoba 0.8.6, and xPerm 1.2.4. The structural tests cover reference and tutorial coverage, metadata, local page links, cell identifiers, parseable inputs, and the contexts of cached outputs. Multiline tutorial cells can contain a list of held input expressions; the input checks accept that standard notebook representation. The tests inspect notebook structure without executing the examples. Evaluate tutorial examples sequentially as described above to check their computational results.

The tutorial was evaluated in a fresh Wolfram Language 15.0.1 session with the same xAct versions: all 337 input cells ran, producing 267 output cells. The only diagnostics were the three intentionally illustrated in the manual (`Validate::inhom`, `Throw::nocatch`, and `SetCMetric::old`). The 13 structural checks pass, and all 41 compiled documentation pages resolve, including the tutorial. Representative compiled excerpts were visually checked. During this headless build, PacletBuild returned Success and created the archive, then its auxiliary front end reported a shutdown fault. A separate session verified all compiled links and rendered the tutorial successfully.

The pages explicitly distinguish declared interfaces from implemented behavior. In this source version, `ZeroDegQ` and `DiffFormQ` have no implementation; the sort order uses `CartanD`; `ExpandHodgeDual` requires a manifold-tagged coframe; and changing `RiemannForm` between two nontrivial connections needs care. These documentation additions do not modify the kernel code.

The file layout, source stylesheets, metadata, and build process follow Wolfram's [Creating Paclets](https://reference.wolfram.com/language/PacletTools/tutorial/CreatingPaclets.html), [Authoring Tech Notes](https://reference.wolfram.com/language/DocumentationTools/tutorial/AuthoringTechNotesUsingDocumentationTools.html), and [PacletDocumentationBuild](https://reference.wolfram.com/language/PacletTools/ref/PacletDocumentationBuild.html) documentation.
