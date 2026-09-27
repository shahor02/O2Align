---
name: code-review
description: Review the current diff, or a PR number/branch/path target, for correctness bugs (plus reuse/simplification/efficiency/decomposition cleanups) at the given effort level (low/medium: fewer, high-confidence findings; high→max: broader coverage, may include uncertain findings); with no level given, it reuses the level you typed last. Pass --comment to post findings as inline PR comments, or --fix to apply the findings to the working tree after the review. O2Align project variant: also flags large, multi-responsibility functions that should be split into logically self-consistent helper methods.
---

# Code review (O2Align)

Follow the standard code-review flow: determine the target (default: uncommitted + unpushed diff on
the current branch; otherwise a given PR number/branch/path), scope the review to changed lines and
their immediately surrounding context, and run at the requested effort level (default: the level last
used in this session, else `medium`).

Review dimensions, in priority order:

1. **Correctness bugs** — logic errors, off-by-one, wrong operator, incorrect ordering
   dependencies (this codebase has several: `Volume::finalise()` level assignment, the
   measurement-leaf-to-root walk in `AlignmentSpec.cxx`), lifetime/ownership issues, and anything
   that would silently produce a wrong alignment result rather than a crash.
2. **Decomposition — monolithic functions.** This is a standing project preference
   ([CLAUDE.md](../../../CLAUDE.md): "Prefer splitting large functions to logically self-consistent
   helper methods."). Flag any changed or touched function that:
   - mixes more than one distinct responsibility in a single body (e.g. parsing + validation +
     computation + I/O all inline), or
   - has grown long enough that a reader must hold several unrelated pieces of state in their head
     at once to follow it, or
   - repeats a block of logic inline that already appears, or could easily be reused, elsewhere in
     the same class/file.
   Prefer this over a purely line-count heuristic: a long but linear, single-purpose loop (especially
   one marked or clearly performance-critical, e.g. per-hit/per-track) is fine as-is per
   CLAUDE.md's efficiency carve-out; a short function that interleaves unrelated concerns is not.
   When flagging, name the specific seams (what the extracted helper(s) would be responsible for),
   not just "this function is too long."
3. **Reuse / simplification / efficiency** — dead code, duplicated logic, unnecessary allocations or
   virtual calls in per-hit/per-track paths, and other cleanups in the standard review recipe.

Do not flag decomposition issues in code you did not touch — this is a review of the diff, not a
drive-by refactor of the file. Do not suggest splitting a function solely because of length if it is
a single coherent loop or an intentionally linear performance-critical path.

Verify every finding (including decomposition ones) before reporting: re-read the function in full
context, confirm the seams you'd extract are real and independently testable/nameable, and rule out
that the "mixed responsibilities" are in fact one cohesive operation that would be harder to follow
split apart.

Report findings with `ReportFindings`, ranked most-severe first. Use category `decomposition` for
monolithic-function findings (separate from `simplification`), with `summary` naming the
responsibilities to split out and `failure_scenario` describing what makes the current shape hard to
verify or maintain (e.g. "a bug in the calibration-derivative branch is easy to miss because it's
interleaved with the rigid-body branch in the same 80-line loop body").

Respect `--comment` (post as inline PR comments) and `--fix` (apply findings to the working tree
after review) exactly as in the standard flow.
