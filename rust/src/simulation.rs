use std::f64::consts::PI;
use rand::Rng;
use rayon::prelude::*;

use crate::testing::example_function;
use crate::types::*;

const FRICTION: f64 = 0.99;
const GRAVITY: [f64; 2] = [0.0, 0.098];
const TARGET_DENSITY: f64 = 2.75;
const PRESSURE_MULTIPLIER: f64 = 0.5;

pub fn set_positions(positions: &mut Vec<[f64; 2]>, window_size: [f64; 2], randomize: bool) {
    if randomize {
        use rand::Rng;
        let mut rng = rand::rng();
        for position in positions.iter_mut() {
            position[0] = rng.random_range(0.0..window_size[0]);
            position[1] = rng.random_range(0.0..window_size[1]);
        }
    } else {
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

pub fn apply_pressure_forces(particle_data: &mut ParticleData, smoothing_radius: f64) {
    let pressure_forces: Vec<[f64; 2]> = (0..particle_data.positions.len())
        .into_par_iter()
        .map(|i| calculate_pressure_force(&particle_data, i, smoothing_radius))
        .collect();
    
        for i in 0..particle_data.velocities.len() {
            if particle_data.densities[i] > 0.0 {
                particle_data.velocities[i][0] += pressure_forces[i][0] / particle_data.densities[i] / particle_data.mass;
                particle_data.velocities[i][1] += pressure_forces[i][1] / particle_data.densities[i] / particle_data.mass;
            }
        }
}

pub fn resolve_collisions(particle_data: &mut ParticleData, window_size: [f64; 2]) {
    for index in 0..particle_data.positions.len() {
        let position = &particle_data.positions[index];
        let velocity = &mut particle_data.velocities[index];

        if position[0] <= 0.0 || position[0] >= window_size[0] {
            velocity[0] = -velocity[0];
        }
        if position[1] <= 0.0 || position[1] >= window_size[1] {
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
    if distance >= radius { return 0.0; }
    //normalized is so that there arent wildly large values messing up the density calculation
    let normalized_radius = radius / 100.0;
    let normalized_distance = distance / 100.0;
    let volume = PI * normalized_radius.powf(4.0) / 6.0;
    let value = normalized_radius - normalized_distance;
    return value * value / volume;
}

pub fn smoothing_kernel_derivative(radius: f64, distance: f64) -> f64 {
    if distance >= radius { return 0.0; }
    let normalized_radius = radius / 100.0;
    let normalized_distance = distance / 100.0;
    let scale = 12.0 / (radius.powf(4.0) * PI);
    return (distance - radius) * scale;
}

pub fn calculate_density(particle_data: &ParticleData, index: usize, smoothing_radius: f64) -> f64 {
    let mut density: f64 = 0.0;
    let reference_position = particle_data.positions[index];

    for j in 0..particle_data.positions.len() {
        if j == index { continue; }
        let other_position = particle_data.positions[j];
        let dx = reference_position[0] - other_position[0];
        let dy = reference_position[1] - other_position[1];
        let distance = (dx * dx + dy * dy).sqrt();
        let influence = smoothing_kernel(smoothing_radius, distance);
        density += particle_data.mass * influence;
    }
    return density;
}

pub fn precalculate_densities(particle_data: &mut ParticleData, smoothing_radius: f64) {
    let densities: Vec<f64> = (0..particle_data.densities.len())
        .into_par_iter()
        .map(|i| calculate_density(&particle_data, i, smoothing_radius))
        .collect();
    
    particle_data.densities = densities;
}

pub fn approximate_value(point: &[f64; 2], particle_data: &ParticleData, smoothing_radius: f64) -> f64 {
    let mut value = 0.0;
    let mass = particle_data.mass;
    let radius_sq = smoothing_radius * smoothing_radius;
    let positions = &particle_data.positions;
    let densities = &particle_data.densities;

    for i in 0..positions.len() {
        let dx = point[0] - positions[i][0];
        let dy = point[1] - positions[i][1];
        let distance_sq = dx * dx + dy * dy;
        
        if distance_sq > radius_sq {
            continue;
        }
        
        let distance = distance_sq.sqrt();
        let influence = smoothing_kernel(smoothing_radius, distance);
        value += example_function(positions[i][0], positions[i][1]) * influence * mass / densities[i];
    }

    return value;
}

pub fn calculate_pressure_force(particle_data: &ParticleData, index: usize, smoothing_radius: f64) -> [f64; 2] {
    let mut net_force = [0.0, 0.0];
    let reference_position = particle_data.positions[index];
    let mut rng = rand::rng();

    for i in 0..particle_data.positions.len() {
        if i == index { continue; }
        let offset: [f64; 2] = [(reference_position[0] - particle_data.positions[i][0]),
                                (reference_position[1] - particle_data.positions[i][1])];
        let distance = (offset[0] * offset[0] + offset[1] * offset[1]).sqrt();
        let direction: [f64; 2] = (|distance: f64| -> [f64; 2] {
            if distance == 0.0 {
                let random_x: f64 = rng.random_range(0.00000000000001..1.0) * 2.0 - 1.0;
                let random_y: f64 = rng.random_range(0.00000000000001..1.0) * 2.0 - 1.0;
                let random_distance = (random_x * random_x + random_y * random_y).sqrt();
                return [random_x / random_distance, random_y / random_distance];
            } else {
                return [offset[0] / distance, offset[1] / distance];
            }
        })(distance);
        let slope = smoothing_kernel_derivative(smoothing_radius, distance);
        let density = particle_data.densities[i];
        net_force[0] += -convert_density_to_pressure(density) * direction[0] * slope * particle_data.mass / density;
        net_force[1] += -convert_density_to_pressure(density) * direction[1] * slope * particle_data.mass / density;
    }

    return net_force;
}

pub fn convert_density_to_pressure(density: f64) -> f64 {
    let density_error = density - TARGET_DENSITY;
    let pressure = density_error * PRESSURE_MULTIPLIER;
    return pressure;
}

pub fn simulation_step(particle_data: &mut ParticleData, smoothing_radius: f64, window_size: [f64; 2]) {
    println!("Simulating step...");
    precalculate_densities(particle_data, smoothing_radius);
    apply_pressure_forces(particle_data, smoothing_radius);
    //apply_gravity(&mut particle_data.velocities);
    apply_friction(&mut particle_data.velocities);
    resolve_collisions(particle_data, window_size);
    update_positions(particle_data);
}