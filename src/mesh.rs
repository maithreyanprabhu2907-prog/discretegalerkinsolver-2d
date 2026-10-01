pub struct Mesh {
    pub x_left: Vec<f64>,
    pub x_right: Vec<f64>,
}

impl Mesh {
    pub fn new(mesh_size: usize) -> Self {
        let mut x_left = vec![0.0; mesh_size];
        let mut x_right = vec![0.0; mesh_size];
        let h = 1.0 / mesh_size as f64;
        for i in 0..mesh_size {
            x_left[i] = i as f64 * h;
            x_right[i] = (1 + i) as f64 * h;
        }
        Mesh { x_left, x_right }
    }

    pub fn n_elements(&self) -> usize {
        self.x_left.len()
    }
    pub fn n_element_width(&self, i: usize) -> f64 {
        self.x_right[i] - self.x_left[i]
    }
}
