use std::f64::consts::PI;

use crate::types::*;

const FRICTION: f64 = 0.99;
const GRAVITY: [f64; 2] = [0.0, 0.098];

pub fn set_positions(positions: &mut Vec<[f64; 2]>, window_size: [f64; 2]) {
    let grid_size = positions.len().isqrt();
    let spacing_x = window_size[0] / (grid_size + 1) as f64;
    let spacing_y = window_size[1] / (grid_size + 1) as f64;
    
    for i in 0..grid_size {
        for j in 0..grid_size {
            let x = spacing_x * (i as f64 + 1.0);
            let y = spacing_y * (j as f64 + 1.0);
            positions[i * grid_size + j] = [x, y];
        }
    }
}

pub fn apply_gravity(velocities: &mut Vec<[f64; 2]>) {
    for velocity in velocities {
        velocity[0] += GRAVITY[0];
        velocity[1] += GRAVITY[1];
    }
}

pub fn apply_friction(velocities: &mut Vec<[f64; 2]>) {
    for velocity in velocities {
        velocity[0] *= FRICTION;
        velocity[1] *= FRICTION;
    }
}

pub fn resolve_collisions(particle_data: &mut ParticleData, window_size: [f64; 2]) {
    for index in 0..particle_data.positions.len() {
        let position = &particle_data.positions[index];
        let velocity = &mut particle_data.velocities[index];

        if position[0] <= 0.0 || position[0] >= window_size[0] {
            velocity[0] = -velocity[0];
        } else if position[1] <= 0.0 || position[1] >= window_size[1] {
            velocity[1] = -velocity[1];
        }
    }
}

pub fn update_positions(particle_data: &mut ParticleData) {
    for index in 0..particle_data.positions.len() {
        particle_data.positions[index][0] += particle_data.velocities[index][0];
        particle_data.positions[index][1] += particle_data.velocities[index][1];
    }
}

pub fn smoothing_kernel(radius: f64, distance: f64) -> f64 {
    let normalized_radius = radius / 100.0;
    let normalized_distance = distance / 100.0;
    let volume = PI * normalized_radius.powf(4.0) / 6.0;
    println!("volume: {}", volume);
    let value = (normalized_radius - normalized_distance).clamp(0.0, normalized_radius);
    println!("Kernel value for distance {}: {}", distance, value);
    return value * value / volume;
}

pub fn calculate_density(particle_data: &ParticleData, index: usize, smoothing_radius: f64) -> f64 {
    let mut density: f64 = 0.0;
    let reference_position = &particle_data.positions[index];

    let mass = 1.0; // Assuming each particle has a mass of 1.0 for simplicity

    for j in 0..particle_data.positions.len() {
        if j == index { continue; }
        let other_position = &particle_data.positions[j];
        let dx = reference_position[0] - other_position[0];
        let dy = reference_position[1] - other_position[1];
        let distance = (dx * dx + dy * dy).sqrt();
        let influence = smoothing_kernel(smoothing_radius, distance);
        density += mass * influence;
        println!("density: {}", density);
    }
    return density;
}