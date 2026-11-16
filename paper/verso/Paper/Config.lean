import VersoPaper.Config
import Paper.References

/-!
# Configuration of the paper

Statement kinds and their numbering, the preamble and the bibliography style.
-/

namespace Paper

open VersoPaper

def config : Config where
  statements := #[
    { name := "theorem", heading := "Theorem", within := some "section" },
    { name := "proposition", heading := "Proposition", counter := some "theorem" },
    { name := "lemma", heading := "Lemma", counter := some "theorem" },
    { name := "corollary", heading := "Corollary", counter := some "theorem" },
    { name := "definition", heading := "Definition", style := "definition", counter := some "theorem" },
    { name := "remark", heading := "Remark", style := "definition", counter := some "theorem" },
    { name := "example", heading := "Example", style := "definition", counter := some "theorem" }]
  notations := #["ob"]
  preamble := include_str "../paper-style.tex"
  jobName := "f4-fsymbols"
  source := "verso/Paper/Text.lean"
  buildCommand := "make -C paper"
  references := references
  bibStyle := "amsplain"
  breakFilePathsFrom := some 36

end Paper
