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

pub fn eval_solution(coeffs: &[f64], mesh: &Mesh, k: usize, x: f64) -> f64 {
    let mut result = 0.0;
    let reference = to_reference(x, mesh.x_left[k], mesh.n_element_width(k));
    for i in 0..=P {
        result += coeffs[(P + 1) * k + i] * legendre(i, reference);
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

pub fn compute_residual(coeffs: &[f64], mesh: &Mesh, res: &mut Vec<f64>) {
    let n = mesh.n_elements();
    for k in 0..n {
        let h = mesh.n_element_width(k);
        for j in 0..=P {
            let mut vol_sum = 0.0;
            for q in 0..3 {
                let xi = GAUSS_POINTS[q];
                let x = to_physical(xi, mesh.x_left[k], h);
                vol_sum +=
                    GAUSS_WEIGHTS[q] * legendre_deriv(j, xi) * eval_solution(coeffs, mesh, k, x);
            }
            let u_right = eval_solution(coeffs, mesh, k, mesh.x_right[k]);

            let u_left = if k == 0 {
                eval_solution(coeffs, mesh, n - 1, mesh.x_right[n - 1])
            } else {
                eval_solution(coeffs, mesh, k - 1, mesh.x_right[k - 1])
            };

            res[(P + 1) * k + j] = (2.0 * j as f64 + 1.0) / h
                * (A * vol_sum - A * (u_right * legendre(j, 1.0) - u_left * legendre(j, -1.0)));
        }
    }
}

pub struct Scratch {
    pub k1: Vec<f64>,
    pub k2: Vec<f64>,
    pub k3: Vec<f64>,
    pub k4: Vec<f64>,
    pub c_temp1: Vec<f64>,
    pub c_temp2: Vec<f64>,
    pub c_temp3: Vec<f64>,
}

impl Scratch {
    pub fn new(size: usize) -> Self {
        Scratch {
            k1: vec![0.0; size],
            k2: vec![0.0; size],
            k3: vec![0.0; size],
            k4: vec![0.0; size],
            c_temp1: vec![0.0; size],
            c_temp2: vec![0.0; size],
            c_temp3: vec![0.0; size],
        }
    }
}

pub fn rk4_step(sol: &mut Solution, mesh: &Mesh, scratch: &mut Scratch, dt: f64) {
    compute_residual(&sol.coefficients, mesh, &mut scratch.k1);
    for (out, (ci, ki)) in scratch
        .c_temp1
        .iter_mut()
        .zip(sol.coefficients.iter().zip(scratch.k1.iter()))
    {
        *out = ci + (dt / 2.0) * ki;
    }

    compute_residual(&scratch.c_temp1, mesh, &mut scratch.k2);
    for (out, (ci, ki)) in scratch
        .c_temp2
        .iter_mut()
        .zip(sol.coefficients.iter().zip(scratch.k2.iter()))
    {
        *out = ci + (dt / 2.0) * ki;
    }

    compute_residual(&scratch.c_temp2, mesh, &mut scratch.k3);
    for (out, (ci, ki)) in scratch
        .c_temp3
        .iter_mut()
        .zip(sol.coefficients.iter().zip(scratch.k3.iter()))
    {
        *out = ci + dt * ki;
    }

    compute_residual(&scratch.c_temp3, mesh, &mut scratch.k4);

    for i in 0..sol.coefficients.len() {
        sol.coefficients[i] += (dt / 6.0)
            * (scratch.k1[i] + 2.0 * scratch.k2[i] + 2.0 * scratch.k3[i] + scratch.k4[i]);
    }
}
