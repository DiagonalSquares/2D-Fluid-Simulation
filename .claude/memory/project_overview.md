---
name: project_overview
description: 2D fluid simulation using SPH (Smoothed Particle Hydrodynamics) with Rust backend and web frontend
type: project
---

## What It Is
A 2D particle-based fluid simulation using SPH (Smoothed Particle Hydrodynamics). Particles interact via density/pressure forces, with smoothing kernels governing how nearby particles influence each other.

## Tech Stack
- **Rust** — core simulation logic, with two build targets:
  - Native desktop app via Piston (windowing/rendering), gated behind `native` feature flag
  - WASM library via wasm-bindgen (for web frontend)
- **TypeScript/Web** — frontend in `web/` (Vite + TS project with `main.ts` and `sim.ts`)
- **Cargo workspace** — root `Cargo.toml` is a workspace with `rust/` as the sole member

## Key Simulation Concepts
- **SPH kernels**: `smoothing_kernel` and `smoothing_kernel_derivative` in `sph.rs`
- **Density**: each particle's density is calculated from nearby particles using the smoothing kernel
- **Pressure**: density is converted to pressure via a linear equation of state (`TARGET_DENSITY`, `PRESSURE_MULTIPLIER`)
- **Simulation loop** (in `simulation.rs`): predict positions -> precalculate densities -> apply pressure forces -> apply friction -> resolve collisions -> update positions

## Constants
- `PARTICLE_COUNT`: 300, `DENSITY_RADIUS`: 100.0 (in `main.rs`)
- `FRICTION`: 0.99, `GRAVITY`: [0.0, 0.098], `TARGET_DENSITY`: 2.0, `PRESSURE_MULTIPLIER`: 0.5 (in `physics.rs`)

## Active Branch Work
- `ai-grid-attempt` branch has work-in-progress grid-based spatial optimization (not on main yet)
