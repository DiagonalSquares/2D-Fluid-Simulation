use piston_window::{Texture, TextureSettings};
use image::{ImageBuffer, Rgba};

mod types;
mod testing;
mod simulation;

use crate::types::*;
use crate::testing::*;
use crate::simulation::*;

const PARTICLE_COUNT: usize = 300;
const RANDOM_INDEX: usize = PARTICLE_COUNT / 2;
const DENSITY_RADIUS: f64 = 100.0;

fn main() {
    let window_size = [1200.0, 800.0];

    let mut window: piston_window::PistonWindow =
        piston_window::WindowSettings::new("Simulation", window_size)
            .exit_on_esc(true)
            .build()
            .unwrap();
    
    let mut particles = ParticleData {
        mass: 1.0,
        positions: vec![[0.0, 0.0]; PARTICLE_COUNT],
        predicted_positions: vec![[0.0, 0.0]; PARTICLE_COUNT],
        velocities: vec![[0.0, 0.0]; PARTICLE_COUNT],
        densities: vec![0.0; PARTICLE_COUNT],
    };
    //println!("test1");
    set_positions(&mut particles.positions, window_size, false);
    //precalculate_densities(&mut particles, DENSITY_RADIUS);
    //calculate_density(&particles, RANDOM_INDEX, DENSITY_RADIUS);
    //println!("test2");

    let mut buffer: ImageBuffer<Rgba<u8>, Vec<u8>> = ImageBuffer::new(window_size[0] as u32, window_size[1] as u32);
    //draw_approximation(&mut buffer, &particles, DENSITY_RADIUS);
    //draw_function(&mut buffer);
    
    let mut texture_context = window.create_texture_context();
    let texture = Texture::from_image(
        &mut texture_context,
        &buffer,
        &TextureSettings::new(),
    ).unwrap();
    
    while let Some(event) = window.next() {
        simulation_step(&mut particles, DENSITY_RADIUS, window_size);

        window.draw_2d(&event, |c, g, _device| {
            piston_window::clear([0.1, 0.1, 0.3, 1.0], g);
            //piston_window::image(&texture, c.transform, g);
            render_particles(&particles, 5.0, c, g);
            //render_particles_from_function(&particles, 5.0, c, g);
            //visualize_density(particles.positions[RANDOM_INDEX], c, g);
        });
    }
}

fn render_particles<G: piston_window::Graphics>(
    particle_data: &ParticleData,
    particle_size: f64,
    c: piston_window::Context,
    g: &mut G,
) {
    let size = particle_size;
    for position in &particle_data.positions {
        let x = position[0];
        let y = position[1];

        let rectangle = piston_window::rectangle::square(x - size, y - size, size * 2.0);
        piston_window::ellipse([1.0, 1.0, 1.0, 1.0], rectangle, c.transform, g);
    }
}

