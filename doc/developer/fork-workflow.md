# Fork workflow

This site documents [`kdoble48/moris`](https://github.com/kdoble48/moris), a fork of
[`kkmaute/moris`](https://github.com/kkmaute/moris), the University of Colorado repository.

## Branches

| Branch | Role |
|---|---|
| `develop` | integration branch; every pull request targets it, and this site is built from it |
| `main` | mirrors the last upstream synchronization; not a development target |
| feature and fix branches | branched from `develop`, merged back through pull requests |

The fork is developed independently. Changes flow from upstream into the fork through the
`sync-develop` workflow; contributions back to upstream are a separate, deliberate decision.

## Pull requests

- Branch from the current `origin/develop`.
- Keep a pull request to one topic.
- Build both trees and run the affected tests (see {doc}`testing`).
- Build-system changes state their measured effect on build time or binary size.
- Cite the bug-ledger ID (below) in the commits of a bug fix.

## The bug ledger

Defects found in MORIS are recorded in `doc/internal/BUG_LEDGER.md` as numbered entries
(`MB-001`, `MB-002`, …). Each entry records:

- **Status** — `open`, `fix-written`, `fixed`, or `merged`.
- **Symptom** — what the user sees, with the exact log signature.
- **Reproducer** — a deck and the commands that show the defect.
- **Code** — the offending lines (`file:line`) for an open entry; the diff of the fix for a fixed
  one.
- **Workaround** — what a deck can do until the fix is merged.

Before diagnosing a crash, a NaN, a frozen run, or a wrong gradient, search the ledger for the
symptom: it may already have a reproducer, a workaround, or a fix branch. When you confirm a new
defect, add the next `MB-NNN` entry immediately, even if the fix comes later. A fix adds a
regression test and reruns the original reproducer before the status advances.

Missing *features* are not ledger entries unless they fail silently.
