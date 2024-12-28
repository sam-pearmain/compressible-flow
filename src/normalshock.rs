#![allow(dead_code)]
use crate::{isentropic::{ valid_specific_heat_ratio, IsentropicFlowError }, numerics::newton_raphson};


pub enum Input {
    UpstreamMachNumber(f64),
    DownstreamMachNumber(f64),
    TemperatureRatio(f64),
    PressureRatio(f64),
    DensityRatio(f64),
    StagnationPressureRatio(f64),
}

pub enum Output {
    UpstreamMachNumber,
    DownstreamMachNumber,
    TemperatureRatio,
    PressureRatio,
    DensityRatio,
    StagnationPressureRatio,
}

struct NormalShock {
    upstream_mach_number: f64,          // M1
    downstream_mach_number: f64,        // M2
    temperature_ratio: f64,             // T2 / T1 (static temperature ratio)
    pressure_ratio: f64,                // p2 / p1
    density_ratio: f64,                 // ρ2 / ρ1
    stagnation_temperature_ratio: f64,  // T02 / T01
    stagnation_pressure_ratio: f64,     // p02 / p01
}

// impl NormalShock {
//     pub fn new(input: Input, output: Output, specific_heat_ratio: Option<f64>) -> Result<NormalShock, IsentropicFlowError> {
//         let specific_heat_ratio = specific_heat_ratio.unwrap_or(1.4);
//         if !valid_specific_heat_ratio(specific_heat_ratio) {
//             return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
//         }

//         Ok(())
//     }
// }

pub fn calculate() {

}

pub fn calc_downstream_mach_from_upstream_mach(upstream_mach: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let downstream_mach_squared: f64 = 
        ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) + 2.0) / 
        (2.0 * specific_heat_ratio * upstream_mach.powi(2) - (specific_heat_ratio - 1.0));
    Ok(downstream_mach_squared.sqrt())
}

pub fn calc_pressure_ratio_from_upstream_mach(upstream_mach: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let pressure_ratio: f64 = 
        (2.0 * specific_heat_ratio * upstream_mach.powi(2) - (specific_heat_ratio - 1.0)) /
        (specific_heat_ratio + 1.0);
    Ok(pressure_ratio)
}

pub fn calc_temperature_ratio_from_upstream_mach(upstream_mach: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let temperature_ratio: f64 = 
        (2.0 * specific_heat_ratio * upstream_mach.powi(2) - (specific_heat_ratio - 1.0)) *
        ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) + 2.0) / 
        ((specific_heat_ratio + 1.0).powi(2) * upstream_mach.powi(2));
    Ok(temperature_ratio)
}

pub fn calc_density_ratio_from_upstream_mach(upstream_mach: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let density_ratio: f64 = 
        ((specific_heat_ratio + 1.0) * upstream_mach.powi(2)) / 
        ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) + 2.0);
    Ok(density_ratio)
}

pub fn calc_stagnation_pressure_ratio_from_upstream_mach(upstream_mach: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let stagnation_pressure_ratio: f64 = 
        calc_density_ratio_from_upstream_mach(upstream_mach, specific_heat_ratio)?.powf(specific_heat_ratio / (specific_heat_ratio - 1.0)) *
        (calc_pressure_ratio_from_upstream_mach(upstream_mach, specific_heat_ratio)? / 1.0).powf(1.0 / (specific_heat_ratio - 1.0));
    Ok(stagnation_pressure_ratio)
}

pub fn calc_upstream_mach_from_pressure_ratio(pressure_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let upstream_mach: f64 = ((pressure_ratio - 1.0) * (specific_heat_ratio + 1.0) / (2.0 * specific_heat_ratio) + 1.0).sqrt();
    Ok(upstream_mach)
}

pub fn calc_upstream_mach_from_temperature_ratio(temperature_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    // solves the quadratic equation in the form aM^2 + bM + c = 0
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let a: f64 = 2.0 * specific_heat_ratio * (specific_heat_ratio - 1.0);
    let b: f64 = 4.0 * specific_heat_ratio 
        - (specific_heat_ratio - 1.0) * (specific_heat_ratio - 1.0) 
        - temperature_ratio * (specific_heat_ratio + 1.0) * (specific_heat_ratio + 1.0);
    let c: f64 = -2.0 * (specific_heat_ratio - 1.0);
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return Err(IsentropicFlowError::WhatTheFuck);
    }
    let upstream_mach: f64 = (-b + discriminant.sqrt()) / (2.0 * a).sqrt();
    Ok(upstream_mach)
}

pub fn calc_upstream_mach_from_density_ratio(density_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let upstream_mach: f64 = (2.0 * density_ratio / (specific_heat_ratio + 1.0 - density_ratio * (specific_heat_ratio - 1.0))).sqrt();
    Ok(upstream_mach)
}

pub fn calc_upstream_mach_from_stagnation_pressure_ratio(stagnation_pressure_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    // uses newton raphson to solve the equation f(M1) = p02/p01
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    if stagnation_pressure_ratio <= 0.0 || stagnation_pressure_ratio > 1.0 {
        return Err(IsentropicFlowError::InvalidPressureRatio);
    }
    
    let f  = |upstream_mach: f64| {
        let alpha: f64 = (specific_heat_ratio + 1.0) * upstream_mach.powi(2)
            / ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) + 2.0);
        let beta: f64 = (specific_heat_ratio + 1.0)
            / (2.0 * specific_heat_ratio * upstream_mach.powi(2) - (specific_heat_ratio - 1.0));
        alpha.powf(specific_heat_ratio / (specific_heat_ratio - 1.0))
            * beta.powf(1.0 / (specific_heat_ratio - 1.0))
            - stagnation_pressure_ratio
    };

    let df = |upstream_mach: f64| {
        let alpha = (specific_heat_ratio + 1.0) * upstream_mach.powi(2)
            / ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) + 2.0);
        let beta = (specific_heat_ratio + 1.0)
            / (2.0 * specific_heat_ratio * upstream_mach.powi(2) - (specific_heat_ratio - 1.0));
        let alpha_derivative = (2.0 / upstream_mach
            - 2.0 * upstream_mach * (specific_heat_ratio - 1.0)
                / ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) + 2.0))
            * alpha;
        let beta_derivative = -4.0 * specific_heat_ratio * upstream_mach * beta
            / (2.0 * specific_heat_ratio * upstream_mach.powi(2)
                - (specific_heat_ratio - 1.0));
        (specific_heat_ratio / (specific_heat_ratio - 1.0))
            * alpha.powf(1.0 / (specific_heat_ratio - 1.0))
            * alpha_derivative
            * beta.powf(1.0 / (specific_heat_ratio - 1.0))
            + alpha.powf(specific_heat_ratio / (specific_heat_ratio - 1.0))
                / (specific_heat_ratio - 1.0)
                * beta.powf((2.0 - specific_heat_ratio) / (specific_heat_ratio - 1.0))
                * beta_derivative
    };

    let upstream_mach = newton_raphson(&f, &df, 2.0, None, None);
    Ok(upstream_mach)
}