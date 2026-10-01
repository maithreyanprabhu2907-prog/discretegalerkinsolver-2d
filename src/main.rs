const P: usize = 2;
const CFL: f64 = 0.5;
const A: f64 = 1.0;
const GAUSS_POINTS: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
const GAUSS_WEIGHTS: [f64; 3] = [0.5555555555555555, 0.8888888888888888, 0.555555555555555555];

struct Mesh {
    x_left: Vec<f64>,
    x_right: Vec<f64>,
}

impl Mesh {
    fn new(mesh_size: usize) -> Self {
        let mut x_left = vec![0.0; mesh_size];
        let mut x_right = vec![0.0; mesh_size];
        let h = 1.0 / mesh_size as f64;
        for i in 0..mesh_size {
            x_left[i] = i as f64 * h;
            x_right[i] = (1 + i) as f64 * h;
        }
        Mesh { x_left, x_right }
    }

    fn n_elements(&self) -> usize {
        self.x_left.len()
    }
    fn n_element_width(&self, i: usize) -> f64 {
        self.x_right[i] - self.x_left[i]
    }
}

struct Solution {
    coefficients: Vec<f64>,
    n_dof: usize,
}
//Solution is where dL and dR are computed
impl Solution {
    fn new(mesh_size: usize) -> Self {
        let coefficients = vec![0.0; mesh_size * (P + 1)];
        Self {
            coefficients,
            n_dof: (P + 1),
        }
    }
}

struct TimeIntegrator {
    t: f64,
    t_end: f64,
    dt: f64,
}

impl TimeIntegrator {
    fn new(t: f64, t_end: f64, h: f64) -> Self {
        TimeIntegrator {
            t,
            t_end,
            dt: Self::cfl_compute(h),
        }
    }

    fn cfl_compute(h: f64) -> f64 {
        CFL * h / (A * (2 * P + 1) as f64)
    }
}

fn legendre(p: usize, xi: f64) -> f64 {
    match p {
        0 => 1.0,
        1 => xi,
        2 => 0.5 * (3.0 * xi * xi - 1.0),
        _ => todo!(),
    }
}

fn legendre_deriv(p: usize, xi: f64) -> f64 {
    match p {
        0 => 0.0,
        1 => 1.0,
        2 => 3.0 * xi,
        _ => todo!(),
    }
}

fn to_reference(x: f64, x_left: f64, h: f64) -> f64 {
    (2.0 / h) * (x - x_left) - 1.0
}

fn eval_solution(sol: &Solution, mesh: &Mesh, k: usize, x: f64) -> f64 {
    let mut result = 0.0;
    let reference = to_reference(x, mesh.x_left[k], mesh.n_element_width(k));
    for i in 0..=P {
        result += sol.coefficients[(P + 1) * k + i] * legendre(i, reference);
    }
    result
}

fn to_physical(xi: f64, x_left: f64, h: f64) -> f64 {
    x_left + (xi + 1.0) * h / 2.0
}

fn sin_ic(x: f64) -> f64 {
    (std::f64::consts::PI * x).sin()
}

fn project_initial_condition(sol: &mut Solution, mesh: &Mesh, u0: fn(f64) -> f64) {
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

fn compute_residual(sol: &Solution, mesh: &Mesh, res: &mut Vec<f64>) {
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

fn main() {
    let mesh_size = 10;
    let mesh = Mesh::new(mesh_size);
    let mut sol = Solution::new(mesh_size);
    let time_integrator = TimeIntegrator::new(0.0, 1.0, mesh.n_element_width(0));
    project_initial_condition(&mut sol, &mesh, sin_ic);

    for k in 0..mesh_size {
        let x_mid = (mesh.x_left[k] + mesh.x_right[k]) / 2.0;
        let u_dg = eval_solution(&sol, &mesh, k, x_mid);
        let u_exact = sin_ic(x_mid);
        println!("k={k}  u_dg={u_dg:.6}  exact={u_exact:.6}");
    }
}
