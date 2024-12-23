pub fn newton_raphson(
    f: &impl Fn(f64) -> f64,
    df: &impl Fn(f64) -> f64,
    x0: f64,
    tolerance: Option<f64>,
) -> f64 {
    // default tolerance 1e-9 unless otherwise given
    let tolerance = tolerance.unwrap_or(1e-9);

    // initialize current root estimate
    let mut x_curr = x0;

    // declare next root estimate and function evaluations
    let mut x_next: f64;
    let mut f_next: f64;
    let mut df_next: f64;

    // evaluate the function and its derivative at the initial guess
    let mut f_curr = f(x_curr);
    let mut df_curr = df(x_curr);

    // iterative solution
    for _ in 0..1000 {
        if df_curr.abs() < 1e-12 {
            panic!("derivative is too small, newton-raphson may fail")
        }
        x_next = x_curr - (f_curr / df_curr);

        // evaluate the function and its derivative at the updated root estimate
        f_next = f(x_next);
        df_next = df(x_next);

        // solver termination on convergence criteria
        if (x_next - x_curr).abs() <= tolerance {
            x_curr = x_next;
            break;
        }

        // store updated values for next iteration
        x_curr = x_next;
        f_curr = f_next;
        df_curr = df_next;
    }

    x_curr
}