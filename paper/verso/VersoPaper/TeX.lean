import VersoPaper.Genre
import VersoPaper.Config

/-!
# TeX rendering

Renders a `Paper` document to LaTeX: the preamble template of the paper, statement environments
generated from its statement kinds, and the body. Numbers of sections, statements, equations and
tables are computed here as well, so that every cross-reference is checked before TeX runs and the
numbering can be compared with the one TeX produces.
-/

open Verso Doc

namespace VersoPaper

/-! ## Numbering -/

inductive LabelKind where
  | section | statement | equation | table
  deriving Repr, BEq, Inhabited

structure LabelInfo where
  kind : LabelKind
  number : String
  deriving Repr, Inhabited

structure NumState where
  sections : Array Nat := #[0, 0, 0]
  statements : Std.HashMap String Nat := {}
  equations : Nat := 0
  tables : Nat := 0
  labels : Std.HashMap String LabelInfo := {}
  order : Array String := #[]

abbrev NumM := ReaderT Config (StateT NumState (Except String))

def sectionNames : Array String := #["section", "subsection", "subsubsection"]

def addLabel (label : String) (info : LabelInfo) : NumM Unit := do
  if (← get).labels.contains label then throw s!"Duplicate label: {label}"
  modify fun s => { s with labels := s.labels.insert label info, order := s.order.push label }

def statementKind (cfg : Config) (name : String) : Except String StatementKind :=
  match cfg.statements.find? (·.name == name) with
  | some k => pure k
  | none => throw s!"Statement kind '{name}' is not declared in the paper's configuration"

/-- The kind whose counter a statement kind uses. -/
def counterOwner (cfg : Config) (k : StatementKind) : Except String StatementKind :=
  match k.counter with
  | none => pure k
  | some c => statementKind cfg c

def currentSectionNumber (level : Nat) : NumM String := do
  let secs := (← get).sections
  return ".".intercalate ((secs.extract 0 level).toList.map toString)

def levelOf (within : String) : Except String Nat :=
  match sectionNames.idxOf? within with
  | some i => pure (i + 1)
  | none => throw s!"Unknown sectioning level: {within}"

partial def numberInline : Inline Paper → NumM Unit
  | .other (.display _ (some label)) _ => do
    modify fun s => { s with equations := s.equations + 1 }
    let n := (← get).equations
    let pre ← match (← read).equationsWithin with
      | none => pure ""
      | some w => do pure ((← currentSectionNumber (← liftExcept (levelOf w))) ++ ".")
    addLabel label { kind := .equation, number := pre ++ toString n }
  | .other _ xs | .emph xs | .bold xs | .link xs _ | .concat xs | .footnote _ xs => xs.forM numberInline
  | _ => pure ()
where
  liftExcept {α} (e : Except String α) : NumM α := match e with
    | .ok a => pure a
    | .error err => throw err

partial def numberBlock : Block Paper → NumM Unit
  | .para xs => xs.forM numberInline
  | .ul items | .ol _ items => items.forM (·.contents.forM numberBlock)
  | .dl items => items.forM fun i => do i.term.forM numberInline; i.desc.forM numberBlock
  | .concat xs | .blockquote xs => xs.forM numberBlock
  | .other (.statement kind label _) xs => do
    let cfg ← read
    let k ← liftExcept (statementKind cfg kind)
    let owner ← liftExcept (counterOwner cfg k)
    if k.numbered then
      let n := ((← get).statements.getD owner.name 0) + 1
      modify fun s => { s with statements := s.statements.insert owner.name n }
      let pre ← match owner.within with
        | none => pure ""
        | some w => do pure ((← currentSectionNumber (← liftExcept (levelOf w))) ++ ".")
      if let some l := label then addLabel l { kind := .statement, number := pre ++ toString n }
    else if label.isSome then throw s!"Unnumbered {kind} cannot have a label"
    xs.forM numberBlock
  | .other (.table spec) xs => do
    modify fun s => { s with tables := s.tables + 1 }
    if let some l := spec.label then addLabel l { kind := .table, number := toString (← get).tables }
    xs.forM numberBlock
  | .other _ xs => xs.forM numberBlock
  | .code _ => pure ()
where
  liftExcept {α} (e : Except String α) : NumM α := match e with
    | .ok a => pure a
    | .error err => throw err

partial def numberPart (p : Part Paper) (level : Nat) : NumM Unit := do
  let md := p.metadata.getD {}
  if level ≥ 1 then
    if md.number then
      if level > sectionNames.size then throw "Too many heading levels"
      modify fun s =>
        let secs := s.sections.set! (level - 1) (s.sections[level - 1]! + 1)
        let secs := (List.range secs.size).foldl (fun a i => if i ≥ level then a.set! i 0 else a) secs
        { s with sections := secs }
      -- Reset the counters numbered within this level.
      let cfg ← read
      for k in cfg.statements do
        if k.within == some sectionNames[level - 1]! then
          modify fun s => { s with statements := s.statements.insert k.name 0 }
      if cfg.equationsWithin == some sectionNames[level - 1]! then
        modify fun s => { s with equations := 0 }
      let n ← currentSectionNumber level
      for l in md.label.toList ++ md.aliases do addLabel l { kind := .section, number := n }
    else if md.label.isSome || !md.aliases.isEmpty then
      throw s!"Unnumbered section '{p.titleString}' cannot have a label"
  p.content.forM numberBlock
  p.subParts.forM (numberPart · (level + 1))

def computeNumbers (cfg : Config) (doc : Part Paper) : Except String NumState := do
  let ((), st) ← ((numberPart doc 0).run cfg).run {}
  return st

/-! ## TeX output -/

structure Context where
  config : Config
  labels : Std.HashMap String LabelInfo
  /-- Content of the `.bbl` file, when the bibliography is to be included in the TeX source. -/
  bbl : Option String := none

abbrev RenderM := ReaderT Context (Except String)

/-- TeX for a character of running text. Non-ASCII characters are written as TeX commands. -/
def charTeX : Char → Except String String
  | '&' => pure "\\&" | '%' => pure "\\%" | '#' => pure "\\#" | '_' => pure "\\_"
  | '$' => pure "\\$" | '{' => pure "\\{" | '}' => pure "\\}"
  | '\\' => pure "\\textbackslash{}" | '^' => pure "\\textasciicircum{}"
  | '~' => pure "\\textasciitilde{}"
  | '\u00A0' => pure "~"         -- no-break space
  | '\u2009' => pure "\\,"       -- thin space
  | '–' => pure "--" | '—' => pure "---"
  | '“' => pure "``" | '”' => pure "''" | '‘' => pure "`" | '’' => pure "'"
  | '…' => pure "\\dots{}"
  | 'ä' => pure "\\\"a" | 'ö' => pure "\\\"o" | 'ü' => pure "\\\"u"
  | 'Ä' => pure "\\\"A" | 'Ö' => pure "\\\"O" | 'Ü' => pure "\\\"U"
  | 'é' => pure "\\'e" | 'è' => pure "\\`e" | 'á' => pure "\\'a" | 'í' => pure "\\'i"
  | 'ó' => pure "\\'o" | 'ç' => pure "\\c{c}" | 'ñ' => pure "\\~n"
  | c =>
    if c.toNat < 128 then pure (String.singleton c)
    else throw s!"No TeX rendering for the character '{c}' (U+{String.ofList (Nat.toDigits 16 c.toNat)})"

def escape (s : String) : RenderM String := do
  let mut out := ""
  for c in s.toList do
    match charTeX c with
    | .ok t => out := out ++ t
    | .error e => throw e
  return out

def command (name text : String) := "\\" ++ name ++ "{" ++ text ++ "}"
def join (xs : Array String) := String.join xs.toList

/-- Text with TeX math between `$` signs, as in metadata strings. -/
def textWithMath (s : String) : RenderM String := do
  let parts := s.splitOn "$"
  if parts.length % 2 == 0 then throw s!"Unbalanced $ in '{s}'"
  let mut out := ""
  for (part, i) in parts.zipIdx do
    out := out ++ (← if i % 2 == 1 then pure ("$" ++ part ++ "$") else escape part)
  return out

def hrefURL (url : String) : String :=
  url.replace "\\" "\\\\" |>.replace "#" "\\#" |>.replace "%" "\\%"

def labelInfo (label : String) : RenderM LabelInfo := do
  match (← read).labels[label]? with
  | some i => pure i
  | none => throw s!"Undefined label: {label}"

partial def plainText : Inline Paper → Option String
  | .text s => some s
  | .code s => some s
  | .concat xs => String.join <$> (xs.toList.mapM plainText)
  | _ => none

/-- The literal text of a role's content: a code element, or plain text. -/
def literal (xs : Array (Inline Paper)) (what : String) : RenderM String :=
  match String.join <$> xs.toList.mapM plainText with
  | some s => pure s
  | none => throw s!"{what} must contain only text or code"

def displayEnv (form : DisplayForm) : String :=
  match form with
  | .single => "equation" | .multline => "multline" | .align => "align" | .gather => "gather"

partial def inlineTeX : Inline Paper → RenderM String
  | .text s => escape s
  | .emph xs => do pure (command "emph" (join (← xs.mapM inlineTeX)))
  | .bold xs => do pure (command "textbf" (join (← xs.mapM inlineTeX)))
  | .code s => do pure (command "texttt" (← escape s))
  | .math .inline s => pure ("$" ++ s ++ "$")
  | .math .display s => pure ("\\[" ++ s ++ "\\]")
  | .linebreak _ => pure "\n"
  | .link xs url => do
    if plainText (.concat xs) == some url then
      pure (command "url" url)
    else
      pure (command "href" (hrefURL url) ++ "{" ++ join (← xs.mapM inlineTeX) ++ "}")
  | .footnote _ xs => do pure (command "footnote" (join (← xs.mapM inlineTeX)))
  | .concat xs => join <$> xs.mapM inlineTeX
  | .image .. => throw "Images are not supported"
  | .other (.cite keys) xs => do
    let refs := (← read).config.references
    for k in keys do
      unless refs.any (·.key == k) do throw s!"Unknown citation: {k}"
    let locator := join (← xs.mapM inlineTeX)
    pure ("\\cite" ++ (if locator.isEmpty then "" else "[" ++ locator ++ "]") ++ "{" ++ ",".intercalate keys ++ "}")
  | .other (.ref label) _ => do
    let info ← labelInfo label
    pure (command (if info.kind == .equation then "eqref" else "ref") label)
  | .other (.display form label) xs => do
    let #[.math .display body] := xs
      | throw "A display role must contain exactly one display formula"
    match form, label with
    | .single, none => pure ("\\[" ++ body ++ "\\]")
    | _, _ =>
      let env := displayEnv form ++ (if label.isNone then "*" else "")
      let lbl := match label with | some l => command "label" l | none => ""
      pure ("\\begin{" ++ env ++ "}" ++ lbl ++ body ++ "\\end{" ++ env ++ "}")
  | .other .file xs => do
    let path ← literal xs "A file name"
    let breakSlashes := (← read).config.breakFilePathsFrom.any (path.length ≥ ·)
    let parts ← (path.splitOn "/").toArray.mapM escape
    -- A break is allowed after a slash that is followed by more of the path.
    let mut out := ""
    for (part, i) in parts.zipIdx do
      if i > 0 then
        out := out ++ "/" ++ (if breakSlashes && part != "" then "\\allowbreak " else "")
      out := out ++ part
    pure (command "texttt" out)
  | .other (.notation name) xs => do
    unless (← read).config.notations.contains name do throw s!"Undeclared notation: {name}"
    -- Code content is literal text (`{ob}`0011``); other content is rendered as inline markup.
    match xs with
    | #[.code s] => pure (command name (← escape s))
    | _ => do pure (command name (join (← xs.mapM inlineTeX)))

/-- The `enumitem` label for a pattern such as `(i)`, `(a)`, `(1)` or `Step 1.`. -/
def listLabel (pattern : String) : RenderM String := do
  let cs := pattern.toList.toArray
  let isWordChar (c : Char) := c.isAlphanum
  let mut out := ""
  let mut found := false
  for h : k in [0:cs.size] do
    let c := cs[k]
    let before := if k == 0 then none else cs[k-1]?
    let after := cs[k+1]?
    let standalone := !(before.any isWordChar) && !(after.any isWordChar)
    let counter := if !standalone || found then none else match c with
      | '1' => some "\\arabic*" | 'i' => some "\\roman*" | 'I' => some "\\Roman*"
      | 'a' => some "\\alph*" | 'A' => some "\\Alph*" | _ => none
    match counter with
    | some t => out := out ++ t; found := true
    | none => out := out ++ (← escape (String.singleton c))
  unless found do throw s!"List label '{pattern}' has no counter (1, i, I, a or A)"
  return out

def columnTeX : Column → String
  | .left => "l" | .center => "c" | .right => "r"
  | .para w => "p{" ++ w ++ "\\linewidth}"

def isListLike : Block Paper → Bool
  | .ul .. | .ol .. | .code .. | .other (.list ..) _ => true
  | _ => false

def isPara : Block Paper → Bool
  | .para .. => true
  | _ => false

/-- Expands `concat` blocks, which Verso uses to group consecutive blocks. -/
partial def flatten (bs : Array (Block Paper)) : Array (Block Paper) :=
  bs.flatMap fun
    | .concat xs => flatten xs
    | b => #[b]

def trimEnd (s : String) : String := s.trimAsciiEnd.toString

mutual

/-- A sequence of blocks. A list or code block belongs to the paragraph around it: no paragraph
break is made before it, after a paragraph, or after it, before a paragraph. -/
partial def blocksTeX (bs : Array (Block Paper)) : RenderM String := do
  let bs := flatten bs
  let mut out := ""
  for h : i in [0:bs.size] do
    let b := bs[i]
    let next := bs[i+1]?
    let s ← blockTeX b
    if isPara b && next.any isListLike then
      out := out ++ trimEnd s ++ "\n"
    else if isListLike b && !(next.any isPara) then
      out := out ++ trimEnd s ++ "\n\n"
    else
      out := out ++ s
  return out

partial def itemsTeX (items : Array (ListItem (Block Paper))) : RenderM String := do
  let mut out := ""
  for item in items do
    out := out ++ "\\item " ++ trimEnd (← blocksTeX item.contents) ++ "\n"
  return out

partial def listTeX (b : Block Paper) (label indent : Option String) : RenderM String := do
  let (env, items) ← match b with
    | .ul items => pure ("itemize", items)
    | .ol 1 items => pure ("enumerate", items)
    | .ol _ _ => throw "Ordered lists must start at 1"
    | _ => throw "A list directive must contain one list"
  let mut opts := #[]
  if let some l := label then
    if env == "itemize" then throw "Unordered lists have no labels"
    opts := opts.push ("label=" ++ (← listLabel l))
  if let some i := indent then opts := opts.push ("leftmargin=" ++ i)
  let optStr := if opts.isEmpty then "" else "[" ++ ",".intercalate opts.toList ++ "]"
  return "\\begin{" ++ env ++ "}" ++ optStr ++ "\n" ++ (← itemsTeX items) ++ "\\end{" ++ env ++ "}\n"

partial def cellTeX (spec : TableSpec) (b : Block Paper) : RenderM (Nat × String) := do
  match b with
  | .para xs => return (1, trimEnd (join (← xs.mapM inlineTeX)))
  | .other (.cell span col) cs =>
    let body ← match flatten cs with
      | #[] => pure ""
      | #[.para xs] => do pure (trimEnd (join (← xs.mapM inlineTeX)))
      | _ => throw s!"A cell must contain at most one paragraph (table {spec.label})"
    let colSpec := match col with
      | some c => columnTeX c
      | none => "l"
    if span == 1 && col.isNone then return (1, body)
    return (span, "\\multicolumn{" ++ toString span ++ "}{" ++ colSpec ++ "}{" ++ body ++ "}")
  | _ => throw s!"Unsupported table cell in table {spec.label}"

partial def rowTeX (spec : TableSpec) (row : ListItem (Block Paper)) : RenderM String := do
  let #[.ul cells] := flatten row.contents
    | throw s!"A table row must be a list of cells (table {spec.label})"
  let mut out := #[]
  let mut width := 0
  for cell in cells do
    let #[c] := flatten cell.contents | throw s!"A table cell must contain one paragraph (table {spec.label})"
    let (span, s) ← cellTeX spec c
    width := width + span
    out := out.push s
  unless width == spec.columns.size do
    throw s!"A row of table {spec.label} spans {width} columns instead of {spec.columns.size}"
  -- Rows that consist of one spanning cell are written without the space before `\\`.
  let sep := if cells.size == 1 && width > 1 then "\\\\" else " \\\\"
  return " & ".intercalate out.toList ++ sep ++ "\n"

partial def tableTeX (spec : TableSpec) (content : Array (Block Paper)) : RenderM String := do
  let cfg := (← read).config
  let mut caption := none
  let mut rows := #[]
  for b in flatten content do
    match b with
    | .other .caption xs => caption := some (trimEnd (← blocksTeX xs))
    | .ul items => rows := rows ++ items
    | _ => throw s!"A table contains a caption and a list of rows (table {spec.label})"
  let mut body := ""
  for h : i in [0:rows.size] do
    if i == spec.header && spec.header > 0 then body := body ++ "\\midrule\n"
    body := body ++ (← rowTeX spec rows[i])
  let mut out := "\\begin{table}[" ++ cfg.tablePlacement ++ "]\n\\centering\n"
  if let some s := spec.size then out := out ++ "\\" ++ s ++ "\n"
  if let some r := spec.rowSpacing then out := out ++ "\\renewcommand{\\arraystretch}{" ++ r ++ "}\n"
  out := out ++ "\\begin{tabular}{" ++ String.join (spec.columns.toList.map columnTeX) ++ "}\n\\toprule\n"
  out := out ++ body ++ "\\bottomrule\n\\end{tabular}\n"
  if let some c := caption then
    out := out ++ cfg.captionSkip ++ "\n" ++ command "caption" c ++ "\n"
  if let some l := spec.label then out := out ++ command "label" l ++ "\n"
  return out ++ "\\end{table}\n\n"

partial def blockTeX : Block Paper → RenderM String
  | .para xs => do pure (trimEnd (join (← xs.mapM inlineTeX)) ++ "\n\n")
  | .concat xs => blocksTeX xs
  | .code s => pure ("\\begin{verbatim}\n" ++ trimEnd s ++ "\n\\end{verbatim}\n")
  | b@(.ul _) | b@(.ol _ _) => listTeX b none none
  | .dl _ => throw "Description lists are not supported"
  | .blockquote xs => do pure ("\\begin{quote}\n" ++ trimEnd (← blocksTeX xs) ++ "\n\\end{quote}\n\n")
  | .other .abstract _ => throw "The abstract must come before the first section"
  | .other (.statement kind label title) xs => do
    let _ ← match statementKind (← read).config kind with
      | .ok k => pure k
      | .error e => throw e
    let opt ← match title with
      | some t => do pure ("[" ++ (← textWithMath t) ++ "]")
      | none => pure ""
    let lbl := match label with | some l => command "label" l | none => ""
    pure ("\\begin{" ++ kind ++ "}" ++ opt ++ lbl ++ "\n" ++ trimEnd (← blocksTeX xs) ++ "\n\\end{" ++ kind ++ "}\n\n")
  | .other (.proof title) xs => do
    let opt ← match title with
      | some t => do pure ("[" ++ (← textWithMath t) ++ "]")
      | none => pure ""
    pure ("\\begin{proof}" ++ opt ++ "\n" ++ trimEnd (← blocksTeX xs) ++ "\n\\end{proof}\n\n")
  | .other (.list label indent) xs => do
    let #[b] := flatten xs | throw "A list directive must contain exactly one list"
    listTeX b label indent
  | .other .paragraph xs => do pure ("\n" ++ (← blocksTeX xs))
  | .other (.table spec) xs => tableTeX spec xs
  | .other .caption _ => throw "A caption must be inside a table"
  | .other (.cell ..) _ => throw "A cell must be inside a table row"
  | .other .bibliography _ => do
    let ctx ← read
    -- The style is declared in both cases: classes such as `amsart` adapt the labels to it.
    let style := "\\bibliographystyle{" ++ ctx.config.bibStyle ++ "}\n"
    match ctx.bbl with
    | some bbl => pure (style ++ trimEnd bbl ++ "\n\n")
    | none => pure (style ++ command "bibliography" ctx.config.jobName ++ "\n\n")

end

partial def partTeX (p : Part Paper) (level : Nat) : RenderM String := do
  let md := p.metadata.getD {}
  let some name := sectionNames[level - 1]? | throw "Too many heading levels"
  let heading := "\\" ++ name ++ (if md.number then "" else "*") ++ "{" ++ join (← p.title.mapM inlineTeX) ++ "}"
  let labels := String.join ((md.label.toList ++ md.aliases).map (command "label"))
  let body ← blocksTeX p.content
  let subs ← p.subParts.mapM (partTeX · (level + 1))
  return heading ++ labels ++ "\n" ++ body ++ join subs

/-- The `\newtheorem` declarations of the statement kinds. -/
def statementDeclarations (cfg : Config) : Except String String := do
  let mut style := "plain"
  let mut out := ""
  for k in cfg.statements do
    if k.style != style then
      out := out ++ "\\theoremstyle{" ++ k.style ++ "}\n"
      style := k.style
    let star := if k.numbered then "" else "*"
    let shared := match k.counter with | some c => "[" ++ c ++ "]" | none => ""
    let within := match k.within with | some w => "[" ++ w ++ "]" | none => ""
    if let some c := k.counter then
      let _ ← statementKind cfg c
    out := out ++ "\\newtheorem" ++ star ++ "{" ++ k.name ++ "}" ++ shared ++ "{" ++ k.heading ++ "}" ++ within ++ "\n"
  return out

def documentTeX (doc : Part Paper) : RenderM String := do
  let cfg := (← read).config
  let md := doc.metadata.getD {}
  let title := join (← doc.title.mapM inlineTeX)
  let shortTitle ← match md.shortTitle with
    | some t => textWithMath t
    | none => pure title
  let authors ← md.authors.mapM fun a => do pure (command "author" (← escape a))
  let date ← match md.date with
    | some d => escape d
    | none => pure ""
  let decls ← match statementDeclarations cfg with
    | .ok d => pure d
    | .error e => throw e
  let preamble := cfg.preamble
    |>.replace "%STATEMENTS%\n" decls
    |>.replace "%TITLE%" title
    |>.replace "%SHORTTITLE%" shortTitle
    |>.replace "%AUTHORS%" ("\n".intercalate authors)
    |>.replace "%DATE%" date
  let (abstract, rest) := (flatten doc.content).partition fun
    | .other .abstract _ => true
    | _ => false
  let abstractTeX ← abstract.mapM fun
    | .other .abstract xs => do pure ("\\begin{abstract}\n" ++ trimEnd (← blocksTeX xs) ++ "\n\\end{abstract}\n")
    | _ => pure ""
  let front := if cfg.abstractBeforeTitle then join abstractTeX ++ "\\maketitle\n\n"
    else "\\maketitle\n\n" ++ join abstractTeX ++ "\n"
  let body ← blocksTeX rest
  let parts ← doc.subParts.mapM (partTeX · 1)
  return s!"% Generated from {cfg.source} by {cfg.buildCommand}. Do not edit.\n" ++
    preamble ++ front ++ body ++ trimEnd (join parts) ++ "\n\n" ++ "\\end{document}\n"

/-! ## BibTeX output -/

def bibTeX (entries : Array BibEntry) : String :=
  String.join <| entries.toList.map fun e =>
    "@" ++ e.kind ++ "{" ++ e.key ++ ",\n" ++
    String.join (e.fields.toList.map fun (k, v) => "  " ++ k ++ " = {" ++ v ++ "},\n") ++ "}\n\n"

/-! ## Labels, for comparison with TeX's -/

def jsonStr (s : String) : String :=
  "\"" ++ (s.replace "\\" "\\\\" |>.replace "\"" "\\\"") ++ "\""

def labelsJson (st : NumState) : String :=
  let kind : LabelKind → String
    | .section => "section" | .statement => "statement" | .equation => "equation" | .table => "table"
  let lines := st.order.toList.map fun l =>
    let i := st.labels[l]!
    "  " ++ jsonStr l ++ ": {\"kind\": " ++ jsonStr (kind i.kind) ++ ", \"number\": " ++ jsonStr i.number ++ "}"
  "{\n" ++ ",\n".intercalate lines ++ "\n}\n"

/-- Checks that every reference resolves, and returns the numbering and the TeX. -/
def render (cfg : Config) (doc : Part Paper) (bbl : Option String) : Except String (NumState × String) := do
  let st ← computeNumbers cfg doc
  let tex ← (documentTeX doc).run { config := cfg, labels := st.labels, bbl }
  return (st, tex)

end VersoPaper
