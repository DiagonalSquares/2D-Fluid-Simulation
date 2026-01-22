use simulation_project::simulation::*;
use simulation_project::types::*;
use image::{ImageBuffer, Rgba};

fn main() {
    let window_size = [800.0, 600.0];

    let mut window: piston_window::PistonWindow =
        piston_window::WindowSettings::new("Simulation", window_size)
            .exit_on_esc(true)
            .build()
            .unwrap();

    let mut particles = ParticleData {
        positions: vec![[0.0, 0.0]; 100],
        velocities: vec![[2.0, 3.0]; 100],
    };
    set_positions(&mut particles.positions, window_size);
    
    while let Some(event) = window.next() {
        //apply_gravity();
        resolve_collisions(&mut particles, window_size);
        apply_friction(&mut particles.velocities);
        update_positions(&mut particles);
        window.draw_2d(&event, |c, g, _device| {
            piston_window::clear([0.1, 0.1, 0.3, 1.0], g);
            render_particles(&particles, 5.0, c, g);
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