# Paper source and rendering

`Paper/Text.lean` is the authoritative source of the paper: prose, mathematics, statements,
proofs, tables, lists and citations, with no layout commands and no raw TeX blocks. Mathematics is
TeX math syntax inside Verso's math markup. `Paper/References.lean` holds the bibliography entries.
`../f4-fsymbols.tex` and `../f4-fsymbols.pdf` are generated: edit the Verso source instead.

Run `make -C paper` from the repository root. This builds the pinned Verso project, renders the
TeX, runs BibTeX on the generated `.bib` to format the cited references, includes the formatted
bibliography in the TeX, compiles it with pdflatex three times from scratch, checks the output
(`check.py`) and copies the TeX and the PDF to `paper/`. Requirements: elan, pdflatex (TeX Live
with `amsart`, `enumitem`, `booktabs`, `mathtools`, `microtype`, `lmodern`), Python 3 and Poppler
(`pdftotext`). A full build takes about 15 seconds once Lake has fetched and built Verso.

`make -C paper submission` also writes `verso/_out/submission/source.zip`: the single TeX file,
which compiles with plain pdflatex runs (no BibTeX), checked by compiling it in an empty
directory. Nothing is submitted. The PDF is reproducible: the build fixes `SOURCE_DATE_EPOCH`.

## Files

- `VersoPaper/`: the genre and the TeX renderer, with nothing specific to this paper.
  - `Genre.lean`: metadata, directives and roles.
  - `Config.lean`: what a paper declares (statement kinds, preamble, references, options).
  - `TeX.lean`: numbering, cross-reference checking, TeX and BibTeX output.
  - `Main.lean`: the command-line driver.
- `Paper/`: this paper.
  - `Text.lean`: the document.
  - `References.lean`: bibliography entries (BibTeX's data model, written to a `.bib` file).
  - `Config.lean`: statement kinds and numbering, bibliography style, options.
  - `Notation.lean`: the paper's text-mode notation role `ob`.
- `paper-style.tex`: the preamble (class, packages, layout, macros) with placeholders for the
  title, authors, date and the generated statement declarations, read when the program runs.
- `build.sh`, `check.py`, `../package.py`, `../Makefile`: build, checks, submission package.
- `lakefile.toml`, `lake-manifest.json`, `lean-toolchain`: Lean `v4.36.0-rc2` and Verso
  `a40afc8c44badcc1e58f060ec70a2aa5cb6243c6`, as in the sheaf-of-rings article.

## Markup

```
# Section title            -- `##` for a subsection; metadata block after the header:
%%%
label := "sec:category"     -- also: aliases := ["sec:other"], number := false
%%%
:::theorem "thm:main"       -- also proposition, lemma, corollary, definition, remark, example,
...                         -- conjecture, question; optional (title := "Name")
:::
::::lemma "lem:x"           -- a title with markup (citations) is a `title` block, first in the statement
:::title
{cite "Rowell"}[Section 3]
:::
...
::::
:::proof (title := "Derivation")   -- title optional
:::
{ref "thm:main"}[]          -- \ref, or \eqref for an equation label
{cite "Rowell" "Gannon"}[Section 3]   -- locator in the brackets, may be empty
$`x^2`   $$`\sum_i x_i`     -- inline and unnumbered display mathematics
{equation "eq:F"}[$$`...`]  -- numbered display; also multline, align, gather (starred without a label)
:::list (label := "(i)") (indent := "2em")   -- around a Verso list; labels (1), (a), (i), Step 1., ...
{file}`code/f4/` {ob}`0011` `inline code`
text[^fn1]   and   [^fn1]: footnote text   -- footnotes
:::listing (size := "footnotesize")          -- around a code block set in a smaller size
:::when "name"  :::unless "name"  {when "name"}[...]  {unless "name"}[...]   -- TeX switches
```

Tables are structured data:

```
::::table "tab:sizes" (columns := "l r p0.4 P0.3") (header := 1) (size := "small") (rowSpacing := "1.05")
:::caption
Sizes of the data.
:::
* * header cell
  * ...
* :::rule
  :::
* * :::cell (span := 3) (column := "p0.9")
    A cell spanning three columns
    :::
::::
```

Rows are list items holding a list of cells; `header` rows are followed by a rule, and a row that
is `:::rule` is a rule between groups of rows; an empty cell is `:::cell` with no content. Columns
are `l`, `c`, `r`, `p<fraction>` (a paragraph column of that fraction of the line width) and
`P<fraction>` (the same, ragged right). A table with the flags `-float -rules` is a centred tabular
without rules inside the running text, with no caption or number. `:::bibliography` places the
bibliography. Non-ASCII typography is written as Unicode and rendered as TeX: no-break space
(U+00A0) as `~`, thin space (U+2009) as `\,`, em space (U+2003) as `\quad`, punctuation space
(U+2008) as `\ `, en and em dashes, curly quotes, `…`, accented letters.

Switches: a paper may declare TeX switches (`\newif\ifname`) in `Config.switches` and set them in
its preamble; `when` and `unless` content renders as `\ifname ... \else ... \fi`, and the renderer
reads the setting so that only typeset content is numbered. This paper declares none.

## Renderer rules

The renderer makes these presentation decisions; everything else follows the markup.

- A list, code block, listing or non-floating table belongs to the paragraph around it: no
  paragraph break is made before it after a paragraph, or after it before a paragraph.
  `:::paragraph` starts a new paragraph after one of them (used once, after the list of
  Theorem 1.1).
- An ellipsis followed by a space is written `\dots\ `, keeping an ordinary interword space; two
  hyphens in typewriter text are written `-{}-`.
- A link whose text is its URL is rendered with `\url`, others with `\href`.
- `{ref}` renders `\eqref` for equation labels and `\ref` otherwise.
- Mathematics in a section heading is written with `\texorpdfstring`, the PDF bookmark showing the
  TeX source without backslashes and braces (two headings here).
- File names are typewriter text; paths of at least `breakFilePathsFrom` characters (36 here) may
  break after a slash (one path in this paper).
- Floating tables are `table` floats (`[ht]`, centred, booktabs rules) with the caption below after
  `\medskip`; only a table with a caption is numbered, as in LaTeX.
- A statement title given as a `title` block is written in braces, `[{...}]`.
- The abstract precedes `\maketitle` (amsart); `\bibliographystyle` is written also when the
  formatted bibliography is included, since amsart adapts its labels to the style.
- The preamble template is read from `paper-style.tex` when the program runs (`preambleFile`), so
  edits to it take effect without rebuilding the Lean code.

Presentation hints that remain in the document: list `indent` (six enumerated lists: the paper
uses different indents for lists with the same label style) and table `size`, `rowSpacing` and
`p` column widths. Itemized lists take their indent from `paper-style.tex`.

## Checks

`check.py` fails on overfull boxes, undefined references, errors or rerun requests in the TeX log;
compares the number of every label computed by the renderer (`_out/pdf/labels.json`) with the one
TeX wrote to the `.aux`; checks that the text ends with the Verso credit (only deferred tables may follow it); and checks that the source
has no TeX commands outside mathematics and no layout commands. The renderer itself fails on an
undefined or duplicate label, an unknown citation key, an undeclared statement kind or notation,
a table row of the wrong width, or an unsupported element.

## Reusing the genre for another paper

`VersoPaper/` contains no paper-specific strings; the same directory, byte for byte, is used by
kim-em/eh-maximal-atlas. For a second paper, copy `VersoPaper/`,
`lakefile.toml`, `lake-manifest.json`, `lean-toolchain`, `Main.lean`, `build.sh` and `check.py`
(changing the job name), and write its own `Paper/` and `paper-style.tex`. Notation roles are
declared with `@[role] meta def name := notationExpander "name"` and listed in
`Config.notations`; the TeX macro of the same name goes in the preamble.
