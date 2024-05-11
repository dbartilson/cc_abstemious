pub mod elements {
    pub mod interpolation {
        pub struct Gp {
            pub coords: [f64; 2],
            pub wt: f64,
        }
        pub static TRIGP: [Gp; 3] = [Gp{coords: [1./6., 1./6.], wt: 1./6.}, 
                                     Gp{coords: [1./6., 2./3.], wt: 1./6.}, 
                                     Gp{coords: [2./3., 1./6.], wt: 1./6.}];
    }   

    pub struct Triangle {

    }
    impl Triangle {
        pub fn shape_functions_at(gp: &interpolation::Gp) -> [f64;3] {
            let xi = gp.coords[0];
            let eta = gp.coords[1];
            let n = [1.0 - xi - eta,
                               xi,
                               eta];
            return n;
        }
    }
    
}