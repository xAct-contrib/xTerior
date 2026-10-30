# xTerior documentation

The editable Wolfram Documentation Tools source notebooks are in:

- `English/ReferencePages/Symbols/`: one page for each of the 39 symbols with an explicit usage message in the two kernel files.
- `English/Guides/xTerior.nb`: the main guide, grouped by topic.

The pages include usage syntax, implementation notes, See Also links, package initialization, self-contained example setup, evaluated examples, and categorization metadata. The documentation link base is `xTerior`.

## Editing

Open a source notebook in Mathematica and use **Palettes > Documentation Tools**. Edit these source notebooks, not the generated notebooks under `build/`. Run each page's example initialization and setup once in a fresh kernel with xAct installed. xTerior's `Diff` is `xAct\`xTerior\`Diff`, distinct from the recent System symbol of the same name.

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

The initial 63 example inputs were evaluated with Wolfram Language 15.0.1, xTensor 1.3.0, xCoba 0.8.6, and xPerm 1.2.4 against the supplied kernel sources. All 102 notebook input cells, including the 39 setup cells, were also evaluated without messages. The ten structural checks cover page completeness, metadata, links, cell identifiers, parseable inputs, and the contexts of cached outputs. All 40 built pages resolve through Documentation Center lookup. The examples are stored in their notebooks, which are the maintained documentation source.

The pages explicitly distinguish declared interfaces from implemented behavior. In this source version, `ZeroDegQ` and `DiffFormQ` have no implementation; the sort order uses `CartanD`; `ExpandHodgeDual` requires a manifold-tagged coframe; and changing `RiemannForm` between two nontrivial connections needs care. These documentation additions do not modify the kernel code.

The file layout, source stylesheets, metadata, and build process follow Wolfram's [Creating Paclets](https://reference.wolfram.com/language/PacletTools/tutorial/CreatingPaclets.html) and [PacletDocumentationBuild](https://reference.wolfram.com/language/PacletTools/ref/PacletDocumentationBuild.html) documentation.
