---
name: project_structure
description: File layout and module responsibilities for the simulation project
type: project
---

## Repo Layout
```
simulation_project/
├── Cargo.toml              # Workspace root, members = ["rust"]
├── Cargo.lock
├── .gitignore              # Ignores /target, /.claude, web/node_modules
├── rust/
│   ├── Cargo.toml          # Package config, edition 2024, cdylib + rlib, "native" feature flag
│   └── src/
│       ├── main.rs         # Native Piston app: window setup, render loop, particle init
│       ├── lib.rs          # Crate root, re-exports all modules; visualization gated behind cfg(feature = "native")
│       ├── types.rs        # ParticleData struct (positions, velocities, densities, mass)
│       ├── simulation.rs   # Simulation loop orchestration (simulation_step)
│       ├── sph.rs          # SPH kernels: smoothing_kernel, smoothing_kernel_derivative, density calculation
│       ├── physics.rs      # Physics: gravity, friction, pressure forces, collision resolution, position updates
│       ├── visualization.rs # Debug/visualization helpers (native-only, behind "native" feature)
│       ├── optimization.rs  # Spatial grid hashing for neighbor lookups (Grid struct, create_grid, create_cell_key, sort_keys)
│       └── api.rs          # WASM API bindings (declared in lib.rs)
└── web/
    ├── package.json
    ├── tsconfig.json
    ├── index.html
    └── src/
        ├── main.ts
        └── sim.ts
```

## Module Dependency Graph
- `main.rs` → uses `simulation_project::*` (the library crate)
- `lib.rs` → re-exports `types`, `sph`, `physics`, `simulation`, `api`; conditionally re-exports `visualization` (native only)
- `simulation.rs` → uses `types`, `sph`, `physics`
- `sph.rs` → uses `types`
- `physics.rs` → uses `types`, `sph`

## Build Targets
- **Native** (default): `cargo build` — includes Piston/OpenGL visualization via `native` feature
- **WASM**: `cargo build --no-default-features --target wasm32-unknown-unknown` — excludes visualization and native deps

## Dependencies (rust/Cargo.toml)
- `wasm-bindgen` — WASM interop
- `rand` — random position initialization
- `rayon` — parallel iteration for density/pressure calculations
- **Optional (native feature):** `piston`, `piston_window`, `piston2d-graphics`, `pistoncore-glutin_window`, `piston2d-opengl_graphics`, `image`
