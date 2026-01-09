use simulation_project::simulation::*;

fn main() {
    let window_size = [800.0, 600.0];

    let mut window: piston_window::PistonWindow =
        piston_window::WindowSettings::new("Simulation", window_size)
            .exit_on_esc(true)
            .build()
            .unwrap();

    let mut simulation = Simulation::new(100, 5.0);
    simulation.set_positions(window_size);
    
    while let Some(event) = window.next() {
        simulation.apply_gravity([0.0, 0.098]);
        simulation.resolve_collisions(window_size);
        simulation.apply_friction();
        simulation.update_positions();
        window.draw_2d(&event, |c, g, _device| {
            piston_window::clear([0.1, 0.1, 0.3, 1.0], g);
            render_particles(&simulation, c, g);
        });
    }
}

fn render_particles<G: piston_window::Graphics>(
    simulation: &Simulation,
    c: piston_window::Context,
    g: &mut G,
) {
    let size = simulation.particle_size;
    for position in &simulation.particle_data.positions {
        let x = position[0];
        let y = position[1];

        let rectangle = piston_window::rectangle::square(x - size, y - size, size * 2.0);
        piston_window::ellipse([1.0, 1.0, 1.0, 1.0], rectangle, c.transform, g);
    }
}