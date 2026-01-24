use simulation_project::simulation::*;
use simulation_project::types::*;

const PARTICLE_COUNT: usize = 100;

fn main() {
    let window_size = [800.0, 600.0];

    let mut window: piston_window::PistonWindow =
        piston_window::WindowSettings::new("Simulation", window_size)
            .exit_on_esc(true)
            .build()
            .unwrap();

    let mut particles = ParticleData {
        positions: vec![[0.0, 0.0]; PARTICLE_COUNT],
        velocities: vec![[0.0, 0.0]; PARTICLE_COUNT],
    };
    set_positions(&mut particles.positions, window_size);

    let index = PARTICLE_COUNT / 2 + 10;
    let radius = 100.0;
    let density = calculate_density(&particles, index, radius);
    println!("Density of particle {}: {}", index, density);
    
    while let Some(event) = window.next() {
        //apply_gravity(&mut particles.velocities);
        resolve_collisions(&mut particles, window_size);
        apply_friction(&mut particles.velocities);
        update_positions(&mut particles);

        window.draw_2d(&event, |c, g, _device| {
            piston_window::clear([0.1, 0.1, 0.3, 1.0], g);
            render_particles(&particles, 5.0, c, g);
            // Draw circle border (outline only) by drawing multiple small circles in a ring
            draw_circle_border(
                particles.positions[index],
                radius,
                2.0,
                [0.0, 0.0, 1.0, 1.0],
                c,
                g,
            );
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

fn draw_circle_border<G: piston_window::Graphics>(
    center: [f64; 2],
    radius: f64,
    width: f64,
    color: [f32; 4],
    c: piston_window::Context,
    g: &mut G,
) {
    // Draw circle border by drawing line segments around the circumference
    let segments = 64;
    let angle_step = std::f64::consts::PI * 2.0 / segments as f64;
    
    for i in 0..segments {
        let angle1 = i as f64 * angle_step;
        let angle2 = ((i + 1) % segments) as f64 * angle_step;
        
        let x1 = center[0] + radius * angle1.cos();
        let y1 = center[1] + radius * angle1.sin();
        let x2 = center[0] + radius * angle2.cos();
        let y2 = center[1] + radius * angle2.sin();
        
        piston_window::line(color, width, [x1, y1, x2, y2], c.transform, g);
    }
}