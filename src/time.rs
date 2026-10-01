use crate::basis::*;

pub struct TimeIntegrator {
    pub t: f64,
    pub t_end: f64,
    pub dt: f64,
}

impl TimeIntegrator {
    pub fn new(t: f64, t_end: f64, h: f64) -> Self {
        TimeIntegrator {
            t,
            t_end,
            dt: Self::cfl_compute(h),
        }
    }

    pub fn cfl_compute(h: f64) -> f64 {
        CFL * h / (A * (2 * P + 1) as f64)
    }
}
