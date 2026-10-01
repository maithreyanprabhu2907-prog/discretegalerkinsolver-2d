pub const P: usize = 2;
pub const CFL: f64 = 0.5;
pub const A: f64 = 1.0;
pub const GAUSS_POINTS: [f64; 3] = [-0.7745966692414834, 0.0, 0.7745966692414834];
pub const GAUSS_WEIGHTS: [f64; 3] = [0.5555555555555555, 0.8888888888888888, 0.555555555555555555];

pub fn to_reference(x: f64, x_left: f64, h: f64) -> f64 {
    (2.0 / h) * (x - x_left) - 1.0
}

pub fn to_physical(xi: f64, x_left: f64, h: f64) -> f64 {
    x_left + (xi + 1.0) * h / 2.0
}

pub fn sin_ic(x: f64) -> f64 {
    (2.0 * std::f64::consts::PI * x).sin()
}
