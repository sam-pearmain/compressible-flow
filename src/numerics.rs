pub fn bisection (
    f: &impl Fn(f64) -> f64,
    x1: f64, // 1st solution bound
    x2: f64, // 2nd solution bound
    tolerance: Option<f64>,
    max_iters: Option<u16>,
) -> f64 {
    // default tolerance 1e-9 unless otherwise given
    let tolerance = tolerance.unwrap_or(1e-9);
    let max_iters: = max_iters.unwrap_or(200);

    // declare midpoint and lower bound and upper bound
    let mut midpoint: f64;
    let mut lowerbound: f64;
    let mut upperbound: f64;

    // initialise lower and upper bound according to given bounds
    if x1 < x2 {
        lowerbound = x1;
        upperbound = x2; 
    } else {
        lowerbound = x2;
        upperbound = x1;
    }

    for _ in 0..max_iters {
        midpoint = (upperbound - lowerbound) / 2.0;
        
        if f(midpoint) < tolerance {
            return midpoint;
        }

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
    let max_iters: = max_iters.unwrap_or(200);
    
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
            break;
        }

        // store updated values for next iteration
        x_current = x_next;
        f_current = f_next;
        df_curr = df_next;
    }

    x_current
}