pub mod basis;
pub mod mesh;
pub mod solution;
pub mod solver;
pub mod time;

pub use basis::sin_ic;
pub use mesh::Mesh;
pub use solution::Solution;
pub use solver::{Scratch, project_initial_condition, rk4_step};
pub use time::TimeIntegrator;

fn main() {}

#[cfg(test)]
mod tests {
    use crate::solver::Scratch;

    use super::*;

    #[test]
    fn test_rk4_full_period() {
        let mesh_size = 100;
        let mesh = Mesh::new(mesh_size);
        let mut sol = Solution::new(mesh_size);
        let mut time = TimeIntegrator::new(0.0, 1.0, mesh.n_element_width(0));
        let mut scratch = Scratch::new(sol.coefficients.len());

        project_initial_condition(&mut sol, &mesh, sin_ic);
        let c0 = sol.coefficients.clone();

        while time.t < time.t_end {
            rk4_step(&mut sol, &mesh, &mut scratch, time.dt);
            time.t += time.dt;
        }

        let max_err = sol
            .coefficients
            .iter()
            .zip(c0.iter())
            .map(|(c_end, c_init)| (c_end - c_init).abs())
            .fold(0.0, f64::max);

        println!("Max coefficient error: {max_err:.6e}");
        assert!(max_err < 0.05, "Error exceeds tolerance: {max_err}");
    }
}
