import VersoPaper.Genre

/-!
# Notation roles of the paper

Text-mode notation used in the prose. Each role is rendered by the TeX macro of the same name,
defined in `paper-style.tex`, which is also available in mathematics.
-/

namespace Paper

open VersoPaper

/-- A simple object, by its Dynkin labels, as in the data files: `{ob}`0011``. -/
@[role] meta def ob := notationExpander "ob"

end Paper
