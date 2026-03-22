---
name: known_issues
description: Remaining organizational issues and potential improvements in the codebase
type: project
---

## Remaining Issues

1. **`api.rs` is empty** — declared in `lib.rs` but has no content yet. Placeholder for WASM API bindings.

2. **Wildcard imports everywhere** — `use crate::module::*` in both `main.rs` and `lib.rs` makes symbol origins hard to trace.

3. **Gravity is commented out** in `simulation_step` — unclear if intentional or forgotten.

4. **Commented-out debug code** scattered in `main.rs` (draw_approximation, draw_function calls).

5. **`web/public/vite.svg`** — leftover Vite scaffold asset, can be deleted if unused.

## Resolved Issues (as of 2026-03-17)
- Simulation.rs split into sph.rs, physics.rs, and simulation.rs modules
- Native graphics deps made optional behind "native" feature flag for WASM compatibility
- testing.rs removed, visualization helpers moved to visualization.rs (native-only)
- .gitignore updated to cover /.claude and web/node_modules
- Duplicate rust/Cargo.lock removed (workspace lock file at root)
