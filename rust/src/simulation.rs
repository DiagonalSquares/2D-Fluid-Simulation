use crate::types::*;

const FRICTION: f64 = 0.99;
const GRAVITY: [f64; 2] = [0.0, 0.098];

pub fn set_positions(positions: &mut Vec<[f64; 2]>, window_size: [f64; 2]) {
    for i in 0..positions.len().isqrt() {
        let x = window_size[0]/positions.len().isqrt() as f64 * i as f64;
        for j in 0..positions.len().isqrt() {
            let y = window_size[1]/positions.len().isqrt() as f64 * j as f64;
            positions[i * 10 + j] = [x + 20.0, y + 20.0];
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