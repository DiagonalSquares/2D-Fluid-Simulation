use crate::types::*;

const FRICTION: f64 = 0.99;

pub struct Simulation {
    pub particle_data: ParticleData,
    pub particle_size: f64,
    pub particle_count: usize,
}

impl Simulation {
    pub fn new(particle_count: usize, particle_size: f64) -> Self {
        Simulation {
            particle_data: ParticleData {
                positions: vec![[0.0; 2]; particle_count],
                velocities: vec![[0.0; 2]; particle_count],
            },
            particle_size,
            particle_count,
        }
    }

    pub fn set_positions(&mut self, window_size: [f64; 2]) {
        for i in 0..self.particle_data.positions.len().isqrt() {
            let x = window_size[0]/self.particle_data.positions.len().isqrt() as f64 * i as f64;
            for j in 0..self.particle_data.positions.len().isqrt() {
                let y = window_size[1]/self.particle_data.positions.len().isqrt() as f64 * j as f64;
                self.particle_data.positions[i * 10 + j] = [x + 20.0, y + 20.0];
            }
        }
    }

    pub fn apply_gravity(&mut self, gravity: [f64; 2]) {
        for velocity in &mut self.particle_data.velocities {
            velocity[0] += gravity[0];
            velocity[1] += gravity[1];
        }
    }

    pub fn apply_friction(&mut self) {
        for velocity in &mut self.particle_data.velocities {
            velocity[0] *= FRICTION;
            velocity[1] *= FRICTION;
        }
    }

    pub fn resolve_collisions(&mut self, window_size: [f64; 2]) {
        for index in 0..self.particle_count {
            let position = &mut self.particle_data.positions[index];
            let velocity = &mut self.particle_data.velocities[index];

            if position[0] <= 0.0 || position[0] >= window_size[0] {
                velocity[0] = -velocity[0];
            } else if position[1] <= 0.0 || position[1] >= window_size[1] {
                velocity[1] = -velocity[1];
            }
        }

    }

    pub fn update_positions(&mut self) {
        for index in 0..self.particle_count {
            self.particle_data.positions[index][0] += self.particle_data.velocities[index][0];
            self.particle_data.positions[index][1] += self.particle_data.velocities[index][1];
        }
    }
}