import VersoPaper.TeX

/-!
# Command-line driver

`paper [--bbl FILE] [--out DIR]` writes `JOB.tex`, `JOB.bib` and `labels.json` to `DIR`
(default `_out`). With `--bbl`, the bibliography is included in the TeX source from that BibTeX
output instead of being read by `\bibliography`, so that the TeX compiles without BibTeX.
-/

open Verso Doc

namespace VersoPaper

def run (cfg : Config) (doc : Part Paper) (args : List String) : IO UInt32 := do
  let rec parse (bbl : Option System.FilePath) (out : System.FilePath) :
      List String → IO (Option System.FilePath × System.FilePath)
    | [] => pure (bbl, out)
    | "--bbl" :: f :: rest => parse (some f) out rest
    | "--out" :: d :: rest => parse bbl d rest
    | a :: _ => throw (IO.userError s!"Unknown argument: {a}")
  let (bblFile, out) ← parse none "_out" args
  let bbl ← bblFile.mapM IO.FS.readFile
  let cfg ← match cfg.preambleFile with
    | some f => do pure { cfg with preamble := ← IO.FS.readFile f }
    | none => pure cfg
  match render cfg doc bbl with
  | .error e =>
    IO.eprintln s!"error: {e}"
    return 1
  | .ok (st, tex) =>
    IO.FS.createDirAll out
    IO.FS.writeFile (out / s!"{cfg.jobName}.tex") tex
    IO.FS.writeFile (out / s!"{cfg.jobName}.bib") (bibTeX cfg.references)
    IO.FS.writeFile (out / "labels.json") (labelsJson st)
    return 0

end VersoPaper
