#[repr(C)]
pub struct ParticleData {
    pub positions: Vec<[f64; 2]>,
    pub velocities: Vec<[f64; 2]>,
}