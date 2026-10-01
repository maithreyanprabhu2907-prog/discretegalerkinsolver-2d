use crate::basis::*;
use crate::mesh::Mesh;
use crate::solution::Solution;

pub fn legendre(p: usize, xi: f64) -> f64 {
    match p {
        0 => 1.0,
        1 => xi,
        2 => 0.5 * (3.0 * xi * xi - 1.0),
        _ => todo!(),
    }
}

pub fn legendre_deriv(p: usize, xi: f64) -> f64 {
    match p {
        0 => 0.0,
        1 => 1.0,
        2 => 3.0 * xi,
        _ => todo!(),
    }
}

pub fn eval_solution(sol: &Solution, mesh: &Mesh, k: usize, x: f64) -> f64 {
    let mut result = 0.0;
    let reference = to_reference(x, mesh.x_left[k], mesh.n_element_width(k));
    for i in 0..=P {
        result += sol.coefficients[(P + 1) * k + i] * legendre(i, reference);
    }
    result
}

pub fn project_initial_condition(sol: &mut Solution, mesh: &Mesh, u0: fn(f64) -> f64) {
    for k in 0..mesh.n_elements() {
        let h = mesh.n_element_width(k);
        for j in 0..=P {
            let mut sum = 0.0;
            for q in 0..3 {
                let xi = GAUSS_POINTS[q];
                let x = to_physical(xi, mesh.x_left[k], h);
                sum += GAUSS_WEIGHTS[q] * u0(x) * legendre(j, xi);
            }
            sol.coefficients[(P + 1) * k + j] = (2 * j + 1) as f64 / 2.0 * sum;
        }
    }
}

pub fn compute_residual(sol: &Solution, mesh: &Mesh, res: &mut Vec<f64>) {
    let n = mesh.n_elements();
    for k in 0..n {
        let h = mesh.n_element_width(k);
        for j in 0..=P {
            let mut vol_sum = 0.0;
            for q in 0..3 {
                let xi = GAUSS_POINTS[q];
                let x = to_physical(xi, mesh.x_left[k], h);
                vol_sum +=
                    GAUSS_WEIGHTS[q] * legendre_deriv(j, xi) * eval_solution(sol, mesh, k, x);
            }
            let u_right = eval_solution(sol, mesh, k, mesh.x_right[k]);

            let u_left = if k == 0 {
                eval_solution(sol, mesh, n - 1, mesh.x_right[n - 1])
            } else {
                eval_solution(sol, mesh, k - 1, mesh.x_right[k - 1])
            };

            res[(P + 1) * k + j] = (2.0 * j as f64 + 1.0) / h
                * (A * vol_sum - A * (u_right * legendre(j, 1.0) - u_left * legendre(j, -1.0)));
        }
    }
}

pub fn rk4_step(sol: &mut Solution, mesh: &Mesh, dt: f64) {
    let size = sol.coefficients.len();
    let mut k1 = vec![0.0; size];
    let mut k2 = vec![0.0; size];
    let mut k3 = vec![0.0; size];
    let mut k4 = vec![0.0; size];

    compute_residual(sol, mesh, &mut k1);
    let mut c_temp1: Vec<f64> = sol
        .coefficients
        .iter()
        .zip(k1.iter())
        .map(|(ci, ki)| ci + (dt / 2.0) * ki)
        .collect();

    let instance1 = Solution {
        coefficients: c_temp1,
        n_dof: P + 1,
    };

    compute_residual(&instance1, mesh, &mut k2);
    let mut c_temp2: Vec<f64> = sol
        .coefficients
        .iter()
        .zip(k2.iter())
        .map(|(ci, ki)| ci + (dt / 2.0) * ki)
        .collect();
    let instance2 = Solution {
        coefficients: c_temp2,
        n_dof: P + 1,
    };

    compute_residual(&instance2, mesh, &mut k3);
    let mut c_temp3: Vec<f64> = sol
        .coefficients
        .iter()
        .zip(k3.iter())
        .map(|(ci, ki)| ci + dt * ki)
        .collect();
    let instance3 = Solution {
        coefficients: c_temp3,
        n_dof: P + 1,
    };
    compute_residual(&instance3, mesh, &mut k4);
    for i in 0..sol.coefficients.len() {
        sol.coefficients[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}
