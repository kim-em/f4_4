/-!
# Paper configuration

What a paper declares to the renderer: its statement kinds and their numbering, its
references, its preamble template and a few presentation choices.
-/

namespace VersoPaper

/-- A theorem-like statement kind. Each kind becomes an `amsthm` environment of the same name. -/
structure StatementKind where
  /-- The directive and environment name, such as `theorem`. -/
  name : String
  /-- The printed heading, such as `Theorem`. -/
  heading : String
  /-- The `amsthm` style: `plain`, `definition` or `remark`. -/
  style : String := "plain"
  /-- Share the counter of this earlier kind. -/
  counter : Option String := none
  /-- Number within sections, such as `section` for numbers 3.1, 3.2, ... -/
  within : Option String := none
  /-- Unnumbered statements. -/
  numbered : Bool := true
  deriving Repr, Inhabited

/-- A bibliography entry, in BibTeX's data model. Field values are BibTeX field text. -/
structure BibEntry where
  kind : String
  key : String
  fields : Array (String × String)
  deriving Repr, Inhabited

structure Config where
  /-- Statement kinds, in declaration order. -/
  statements : Array StatementKind
  /-- The text-mode notation roles the paper declares, by TeX macro name. -/
  notations : Array String := #[]
  /-- TeX switches (`\newif\ifname`) used by `when` and `unless`, by name. The preamble declares
  each one and sets it with `\nametrue` or `\namefalse`; the renderer reads the setting from there,
  so that only the content TeX will typeset is numbered. -/
  switches : Array String := #[]
  /-- The preamble template: everything before the document body. Placeholders:
  `%STATEMENTS%`, `%TITLE%`, `%SHORTTITLE%`, `%AUTHORS%`, `%DATE%`. -/
  preamble : String := ""
  /-- A file to read the preamble template from when the program runs, relative to the working
  directory. Lake does not rebuild a module when a file it reads with `include_str` changes, so
  reading the file at run time makes every edit of it take effect. -/
  preambleFile : Option System.FilePath := none
  /-- Base name of the generated TeX, PDF and BibTeX files. -/
  jobName : String
  /-- Path of the document source, for the header of the generated TeX. -/
  source : String
  /-- The command that regenerates the output, for the header of the generated TeX. -/
  buildCommand : String
  references : Array BibEntry
  /-- The BibTeX style. -/
  bibStyle : String := "plain"
  /-- Whether the abstract precedes `\maketitle`, as `amsart` requires. -/
  abstractBeforeTitle : Bool := true
  /-- Float placement of tables. -/
  tablePlacement : String := "ht"
  /-- Vertical space between a table and its caption. -/
  captionSkip : String := "\\medskip"
  /-- Equations are numbered within this sectioning level, or consecutively. -/
  equationsWithin : Option String := none
  /-- Allow line breaks after the slashes of file paths of at least this many characters. -/
  breakFilePathsFrom : Option Nat := none

end VersoPaper
