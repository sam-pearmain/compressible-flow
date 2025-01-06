pub fn bisection (
    f: &impl Fn(f64) -> f64,
    x1: f64, // 1st solution bound
    x2: f64, // 2nd solution bound
    tolerance: Option<f64>,
    max_iters: Option<u16>,
) -> f64 {
    // default tolerance 1e-9 unless otherwise given
    let tolerance = tolerance.unwrap_or(1e-9);
    let max_iters = max_iters.unwrap_or(200);

    // initialise lower and upper bound according to given bounds
    let (mut lowerbound, mut upperbound) = if x1 < x2 { (x1, x2) } else { (x2, x1) };

    // iterate
    for _ in 0..max_iters {
        let midpoint = (upperbound + lowerbound) / 2.0;
        
        // check convergence
        if f(midpoint).abs() < tolerance || (upperbound - lowerbound) / 2.0 < tolerance {
            return midpoint;
        }

        // update bounds
        if (f(midpoint) * f(lowerbound)) > 0.0 {
            lowerbound = midpoint;
        } else {
            upperbound = midpoint;
        }
    }

    panic!("solution not converged");
}

pub fn newton_raphson(
    f: &impl Fn(f64) -> f64,
    df: &impl Fn(f64) -> f64,
    x_init: f64,
    tolerance: Option<f64>,
    max_iters: Option<u16>,
) -> f64 {
    // default tolerance 1e-9 unless otherwise given
    let tolerance = tolerance.unwrap_or(1e-9);
    let max_iters = max_iters.unwrap_or(200);
    
    // declare next root estimate and function evaluations
    let mut x_current = x_init;
    let mut x_next: f64;
    let mut f_next: f64;
    let mut df_next: f64;

    // evaluate the function and its derivative at the initial guess
    let mut f_current = f(x_current);
    let mut df_curr = df(x_current);

    // iterative solution
    for _ in 0..max_iters {
        if df_curr.abs() < 1e-12 {
            panic!("derivative is too small, newton-raphson may fail")
        }
        x_next = x_current - (f_current / df_curr);

        // evaluate the function and its derivative at the updated root estimate
        f_next = f(x_next);
        df_next = df(x_next);

        // solver termination on convergence criteria
        if (x_next - x_current).abs() <= tolerance {
            x_current = x_next;
            return x_current;
        }

        // store updated values for next iteration
        x_current = x_next;
        f_current = f_next;
        df_curr = df_next;
    }

    panic!("solution not converged");
}

pub fn runge_kutta_fourth_order() {
    
}

pub fn runge_kutta_fourth_order_step(f: &impl Fn(f64, f64, f64) -> f64, r: f64, y1: f64, y2: f64, h: f64) -> (f64, f64) {
    let k1_y1: f64 = h * y2;
    let k1_y2: f64 = h * f(r, y1, y2);

    let k2_y1: f64 = h * (y2 + 0.5 * k1_y2);
    let k2_y2: f64 = h * f(r + 0.5 * h, y1 + 0.5 * k1_y1, y2 + 0.5 * k1_y2);

    let k3_y1: f64 = h * (y2 + 0.5 * k2_y2);
    let k3_y2: f64 = h * f(r + 0.5 * h, y1 + 0.5 * k2_y1, y2 + 0.5 * k2_y2);

    let k4_y1: f64 = h * (y2 + k3_y2);
    let k4_y2: f64 = h * f(r + h, y1 + k3_y1, y2 + k3_y2);

    let new_y1: f64 = y1 + (k1_y1 + 2.0 * k2_y1 + 2.0 * k3_y1 + k4_y1) / 6.0;
    let new_y2: f64 = y2 + (k1_y2 + 2.0 * k2_y2 + 2.0 * k3_y2 + k4_y2) / 6.0;

    (new_y1, new_y2)
}