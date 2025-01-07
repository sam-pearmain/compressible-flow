#![allow(dead_code)]

use crate::isentropic::*; 

pub enum Input {
    UpstreamMach(f64),
    ConeAngle(f64),
    ShockAngle(f64),
    SurfaceMachNumber(f64),
}

pub enum Output {
    SurfaceMachNumber,
    ConeAngle,
    ShockAngle,
    ShockTurnAngle,
    PressureRatio,
    DensityRatio,
    TemperatureRatio,
    StagnationPressureRatio,
    SurfacePressureRatio,
    SurfaceDensityRatio, 
    SurfaceTemperatureRatio,        
    SurfaceStagnationPressureratio, 
}

pub struct SupersonicCone {
    surface_mach_number: f64,               // Mc (mach number at the surface of the cone)
    cone_angle: f64,                        // σ
    shock_angle: f64,                       // β
    shock_turn_angle: f64,                  // δ
    pressure_ratio: f64,                    // p2 / p1
    density_ratio: f64,                     // ρ2 / ρ1
    temperature_ratio: f64,                 // T2 / T1
    stagnation_pressure_ratio: f64,         // p02 / p01
    surface_pressure_ratio: f64,            // pc / p1
    surface_density_ratio: f64,             // ρc / ρ1
    surface_temperature_ratio: f64,         // Tc / T1
    surface_stagnation_pressure_ratio: f64, // p0c / p01 (will equal p02 / p01 since perfectly isentropic compression)
}

impl SupersonicCone {
    pub fn from_mach_and_cone_angle() {
        
    }

    pub fn from_mach_and_shock_angle() {
        
    }

    pub fn from_mach_and_surface_mach() {
        
    }
}

pub fn solve_taylor_maccoll(
    initial_velocity_vector: (f64, f64),
    initial_angle: f64, // typically the shock angle created by the cone
    final_angle: f64,   // the final integration bound
    specific_heat_ratio: f64,
    // a bool to halt the integration once the tangential velocity becomes 0
    steps: Option<i32>,
) -> Result<(Vec<(f64, f64)>, Vec<f64>), IsentropicFlowError> {
    // uses a 4th order runge-kutta method to solve the taylor maccoll equations
    // between the two bounds of initial_angle and final_angle
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }

    // set up step size
    let steps: i32 = steps.unwrap_or(2000);
    let h: f64 = (initial_angle - final_angle) / steps as f64;

    // vectors to store the results
    let mut velocity_components: Vec<(f64, f64)> = Vec::new();
    velocity_components.push(initial_velocity_vector);

    let thetas: Vec<f64> = (0..=steps)
        .map(|i| initial_angle + i as f64 * h)
        .collect();

    // set up starting conditions
    let mut current_radial_velocity: f64 = initial_velocity_vector.0;
    let mut current_tangential_velocity: f64 = initial_velocity_vector.1;
    
    for &theta in &thetas {
        // first runge-kutta constants
        let k1: (f64, f64) = 
            taylor_maccoll(
                (current_radial_velocity, current_tangential_velocity),
                theta,
                specific_heat_ratio,
            )?;
        let k1_radial: f64 = h * k1.0;
        let k1_tangential: f64 = h * k1.1;

        // second runge-kutta constants
        let k2: (f64, f64) =
            taylor_maccoll(
                (current_radial_velocity + (0.5 * k1_radial), current_tangential_velocity + (0.5 * k1_tangential)),
                theta + (0.5 * h),
                specific_heat_ratio,
            )?;
        let k2_radial: f64 = h * k2.0;
        let k2_tangential: f64 = h * k2.1;

        // third runge-kutta constants
        let k3: (f64, f64) =
            taylor_maccoll(
                (current_radial_velocity + (0.5 * k2_radial), current_tangential_velocity + (0.5 * k2_tangential)),
                theta + (0.5 * h),
                specific_heat_ratio,
            )?;
        let k3_radial: f64 = h * k3.0;
        let k3_tangential: f64 = h * k3.1;

        // fourth runge-kutta constants
        let k4: (f64, f64) = 
            taylor_maccoll(
                (current_radial_velocity + k3_radial, current_tangential_velocity + k3_tangential),
                theta + h,
                specific_heat_ratio,
            )?;
        let k4_radial: f64 = h * k4.0;
        let k4_tangential: f64 = h * k4.1;

        // calculate subsequent radial and tangential velocity
        let next_radial_velocity: f64 = 
            current_radial_velocity + (1.0 / 6.0) *
            (k1_radial + 2.0 * k2_radial + 2.0 * k3_radial + k4_radial);
        let next_tangential_velocity: f64 = 
            current_tangential_velocity + (1.0 / 6.0) *
            (k1_tangential + 2.0 * k2_tangential + 2.0 * k3_tangential + k4_tangential);

        // append velocity components to results vector
        velocity_components.push((next_radial_velocity, next_tangential_velocity));
        
        // put break clause here

        // update current values with subsequent values and loop
        current_radial_velocity = next_radial_velocity;
        current_tangential_velocity = next_tangential_velocity;
    }

    Ok((velocity_components, thetas))
}

pub fn taylor_maccoll(velocity_vector: (f64, f64), theta: f64, specific_heat_ratio: f64) -> Result<(f64, f64), IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }

    // extract the radial and tangential velocity components
    let radial_velocity: f64 = velocity_vector.0;
    let tangential_velocity: f64 = velocity_vector.1;

    // calculate the radial velocity derivative
    let radial_velocity_derivative = tangential_velocity;

    // caclulate the tangential velocity derivative
    let denominator: f64 = 
        ((specific_heat_ratio - 1.0) / 2.0) 
        * (1.0 - radial_velocity.powi(2) - tangential_velocity.powi(2)) 
        - tangential_velocity.powi(2);
    let numerator: f64 = 
        tangential_velocity.powi(2) * radial_velocity
        - ((specific_heat_ratio - 1.0) / 2.0)
        * (1.0 - radial_velocity.powi(2) - tangential_velocity.powi(2))
        * (2.0 * radial_velocity + tangential_velocity / theta.tan());

    let tangential_velocity_derivative = numerator / denominator;

    Ok((radial_velocity_derivative, tangential_velocity_derivative))
}