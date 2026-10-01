use crate::basis::P;

pub struct Solution {
    pub coefficients: Vec<f64>,
    pub n_dof: usize,
}
impl Solution {
    pub fn new(mesh_size: usize) -> Self {
        let coefficients = vec![0.0; mesh_size * (P + 1)];
        Self {
            coefficients,
            n_dof: (P + 1),
        }
    }
}
