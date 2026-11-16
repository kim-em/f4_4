import Verso

/-!
# A genre for mathematics papers

Semantic blocks and inlines for a mathematics paper: front matter, sections with labels,
numbered statements and proofs, numbered displays, tables as structured data, lists with
labels, citations and cross-references. Nothing here is specific to one paper: the kinds of
statements, their headings and numbering, and any notation roles are declared by the paper.
-/

open Verso Doc Elab ArgParse
open Lean (quote)

namespace VersoPaper

/-- Metadata of the whole document (root part) and of sections. -/
structure Metadata where
  /-- Authors, root part only. -/
  authors : List String := []
  /-- Date as printed, root part only. -/
  date : Option String := none
  /-- Running head, root part only. TeX math between `$` signs is allowed. -/
  shortTitle : Option String := none
  /-- Cross-reference label of a section. -/
  label : Option String := none
  /-- Further labels of the same section. -/
  aliases : List String := []
  /-- Whether the section is numbered. -/
  number : Bool := true
  deriving Inhabited, Repr, TypeName

/-- The layout of a displayed formula with several lines. -/
inductive DisplayForm where
  /-- One formula. -/
  | single
  /-- One long formula broken over several lines (`\\` separates the lines). -/
  | multline
  /-- Several formulas aligned at `&`. -/
  | align
  /-- Several centred formulas. -/
  | gather
  deriving Repr, Inhabited, BEq

/-- A table column. -/
inductive Column where
  | left | center | right
  /-- A paragraph column whose width is the given fraction of the line width. -/
  | para (width : String)
  deriving Repr, Inhabited, BEq

structure TableSpec where
  label : Option String
  columns : Array Column
  /-- Number of header rows, separated from the body by a rule. -/
  header : Nat := 0
  /-- Font size of the table body (`small`, `footnotesize`, `scriptsize`), if not normal. -/
  size : Option String := none
  /-- Row spacing factor, if not 1. -/
  rowSpacing : Option String := none
  deriving Repr, Inhabited

inductive BlockKind where
  | abstract
  /-- A theorem-like statement. `kind` is one of the statement kinds declared by the paper. -/
  | statement (kind : String) (label : Option String) (title : Option String)
  | proof (title : Option String)
  /-- A list whose items are labelled by `label`, a pattern such as `(i)`, `(a)`, `(1)` or
  `Step 1.`; `indent` is the width reserved for the labels, if not the default. -/
  | list (label : Option String) (indent : Option String)
  /-- Starts a new paragraph where a paragraph would otherwise continue after a list or code. -/
  | paragraph
  | table (spec : TableSpec)
  | caption
  /-- A table cell spanning `span` columns, with its own column format. -/
  | cell (span : Nat) (column : Option Column)
  /-- The bibliography, generated from the paper's references. -/
  | bibliography
  deriving Repr, Inhabited

inductive InlineKind where
  /-- A citation of one or more references; the content is the locator, possibly empty. -/
  | cite (keys : List String)
  /-- A cross-reference to a section, statement, table or equation by label. -/
  | ref (label : String)
  /-- A displayed formula, numbered when it has a label; the content is display math. -/
  | display (form : DisplayForm) (label : Option String)
  /-- A file or directory name. -/
  | file
  /-- A text-mode notation declared by the paper, rendered by the TeX macro of that name. -/
  | notation (name : String)
  deriving Repr, Inhabited

def Paper : Genre where
  PartMetadata := Metadata
  Block := BlockKind
  Inline := InlineKind
  TraverseContext := Unit
  TraverseState := Unit

instance : Repr Paper.Block := inferInstanceAs (Repr BlockKind)
instance : Repr Paper.Inline := inferInstanceAs (Repr InlineKind)
instance : TypeName Paper.PartMetadata := inferInstanceAs (TypeName Metadata)
instance : Inhabited Paper.PartMetadata := inferInstanceAs (Inhabited Metadata)

instance : Lean.Quote Column where
  quote
    | .left => Lean.Syntax.mkCApp ``Column.left #[]
    | .center => Lean.Syntax.mkCApp ``Column.center #[]
    | .right => Lean.Syntax.mkCApp ``Column.right #[]
    | .para w => Lean.Syntax.mkCApp ``Column.para #[quote w]

instance : Lean.Quote DisplayForm where
  quote
    | .single => Lean.Syntax.mkCApp ``DisplayForm.single #[]
    | .multline => Lean.Syntax.mkCApp ``DisplayForm.multline #[]
    | .align => Lean.Syntax.mkCApp ``DisplayForm.align #[]
    | .gather => Lean.Syntax.mkCApp ``DisplayForm.gather #[]

/-- Parses a column: `l`, `c`, `r`, or `p0.4` for a paragraph column of 0.4 line widths. -/
def Column.parse (s : String) : Except String Column :=
  match s with
  | "l" => pure .left
  | "c" => pure .center
  | "r" => pure .right
  | _ =>
    if s.startsWith "p" && (s.drop 1).toString.all (fun c => c.isDigit || c == '.') && s.length > 1 then
      pure (.para (s.drop 1).toString)
    else throw s!"Unknown column '{s}': expected l, c, r or p<fraction>"

def parseColumns (s : String) : Except String (Array Column) :=
  (s.splitOn " " |>.filter (· ≠ "") |>.toArray).mapM Column.parse

/-! ## Arguments -/

structure StatementArgs where
  label : Option String
  title : Option String

instance : FromArgs StatementArgs DocElabM where
  fromArgs := StatementArgs.mk <$> optional (.positional `label .string) <*> .named `title .string true

structure ProofArgs where
  title : Option String

instance : FromArgs ProofArgs DocElabM where
  fromArgs := ProofArgs.mk <$> .named `title .string true

structure ListArgs where
  label : Option String
  indent : Option String

instance : FromArgs ListArgs DocElabM where
  fromArgs := ListArgs.mk <$> .named `label .string true <*> .named `indent .string true

structure TableArgs where
  label : Option String
  columns : String
  header : Nat
  size : Option String
  rowSpacing : Option String

instance : FromArgs TableArgs DocElabM where
  fromArgs := TableArgs.mk <$> optional (.positional `label .string) <*> .named `columns .string false <*>
    .namedD `header .nat 0 <*> .named `size .string true <*> .named `rowSpacing .string true

structure CellArgs where
  span : Nat
  column : Option String

instance : FromArgs CellArgs DocElabM where
  fromArgs := CellArgs.mk <$> .namedD `span .nat 1 <*> .named `column .string true

structure LabelArg where
  label : Option String

instance : FromArgs LabelArg DocElabM where
  fromArgs := LabelArg.mk <$> optional (.positional `label .string)

structure CiteArgs where
  keys : List String

instance : FromArgs CiteArgs DocElabM where
  fromArgs := CiteArgs.mk <$> .many (.positional `key .string)

structure RefArg where
  label : String

instance : FromArgs RefArg DocElabM where
  fromArgs := RefArg.mk <$> .positional `label .string

/-! ## Directives -/

@[directive]
meta def abstract : DirectiveExpanderOf Unit := fun () content => do
  let content ← content.mapM elabBlock
  ``(Block.other (genre := Paper) BlockKind.abstract #[$content,*])

/-- The expander of the directive for a statement kind. -/
meta def statementExpander (kind : String) : DirectiveExpanderOf StatementArgs := fun args content => do
  let content ← content.mapM elabBlock
  ``(Block.other (genre := Paper) (BlockKind.statement $(quote kind) $(quote args.label) $(quote args.title)) #[$content,*])

@[directive]
meta def «theorem» := statementExpander "theorem"
@[directive]
meta def proposition := statementExpander "proposition"
@[directive]
meta def lemma := statementExpander "lemma"
@[directive]
meta def corollary := statementExpander "corollary"
@[directive]
meta def conjecture := statementExpander "conjecture"
@[directive]
meta def definition := statementExpander "definition"
@[directive]
meta def remark := statementExpander "remark"
@[directive]
meta def «example» := statementExpander "example"
@[directive]
meta def question := statementExpander "question"

@[directive]
meta def proof : DirectiveExpanderOf ProofArgs := fun args content => do
  let content ← content.mapM elabBlock
  ``(Block.other (genre := Paper) (BlockKind.proof $(quote args.title)) #[$content,*])

@[directive]
meta def list : DirectiveExpanderOf ListArgs := fun args content => do
  let content ← content.mapM elabBlock
  ``(Block.other (genre := Paper) (BlockKind.list $(quote args.label) $(quote args.indent)) #[$content,*])

@[directive]
meta def paragraph : DirectiveExpanderOf Unit := fun () content => do
  let content ← content.mapM elabBlock
  ``(Block.other (genre := Paper) BlockKind.paragraph #[$content,*])

@[directive]
meta def table : DirectiveExpanderOf TableArgs := fun args content => do
  let columns ← match parseColumns args.columns with
    | .ok cs => pure cs
    | .error e => throwError e
  let content ← content.mapM elabBlock
  ``(Block.other (genre := Paper)
      (BlockKind.table (TableSpec.mk $(quote args.label) $(quote columns) $(quote args.header) $(quote args.size) $(quote args.rowSpacing)))
      #[$content,*])

@[directive]
meta def caption : DirectiveExpanderOf Unit := fun () content => do
  let content ← content.mapM elabBlock
  ``(Block.other (genre := Paper) BlockKind.caption #[$content,*])

@[directive]
meta def cell : DirectiveExpanderOf CellArgs := fun args content => do
  let column ← match args.column.map Column.parse with
    | none => pure none
    | some (.ok c) => pure (some c)
    | some (.error e) => throwError e
  let content ← content.mapM elabBlock
  ``(Block.other (genre := Paper) (BlockKind.cell $(quote args.span) $(quote column)) #[$content,*])

@[directive]
meta def bibliography : DirectiveExpanderOf Unit := fun () content => do
  unless content.isEmpty do throwError "The bibliography is generated from the references; it takes no content"
  ``(Block.other (genre := Paper) BlockKind.bibliography #[])

/-! ## Roles -/

@[role]
meta def cite : RoleExpanderOf CiteArgs := fun args content => do
  if args.keys.isEmpty then throwError "cite needs at least one key"
  let content ← content.mapM elabInline
  ``(Inline.other (genre := Paper) (InlineKind.cite $(quote args.keys)) #[$content,*])

@[role]
meta def ref : RoleExpanderOf RefArg := fun arg content => do
  unless content.isEmpty do throwError "ref takes no content: the number is generated"
  ``(Inline.other (genre := Paper) (InlineKind.ref $(quote arg.label)) #[])

/-- The expander of a role for a displayed formula of the given form. -/
meta def displayExpander (form : DisplayForm) : RoleExpanderOf LabelArg := fun arg content => do
  let content ← content.mapM elabInline
  ``(Inline.other (genre := Paper) (InlineKind.display $(quote form) $(quote arg.label)) #[$content,*])

@[role]
meta def equation := displayExpander .single
@[role]
meta def multline := displayExpander .multline
@[role]
meta def align := displayExpander .align
@[role]
meta def gather := displayExpander .gather

@[role]
meta def file : RoleExpanderOf Unit := fun () content => do
  let content ← content.mapM elabInline
  ``(Inline.other (genre := Paper) InlineKind.file #[$content,*])

/-- The expander of a paper-specific notation role, rendered by the TeX macro `name`. -/
meta def notationExpander (name : String) : RoleExpanderOf Unit := fun () content => do
  let content ← content.mapM elabInline
  ``(Inline.other (genre := Paper) (InlineKind.notation $(quote name)) #[$content,*])

end VersoPaper
