use std::f64::consts::PI;
use crate::numerics;

#[derive(Debug)]
pub enum IsentropicFlowError {
    InvalidMachNumber,
    InvalidMachAngle,
    InvalidTemperatureRatio,
    InvalidPressureRatio,
    InvalidDensityRatio,
    InvalidPrandtlMeyerAngle,
    InvalidSpecificHeatRatio,
}

pub fn calc_mach_from_mach_angle(mach_angle: f64) -> Result<f64, IsentropicFlowError> {
    if mach_angle < 0.0 || mach_angle > PI / 2.0 {
        // check valid mach angle in radians
        return Err(IsentropicFlowError::InvalidMachAngle)
    }
    let mach_number: f64 = 1.0 / mach_angle.sin();
    Ok(mach_number)
}

pub fn calc_mach_from_temperature_ratio(temperature_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    if temperature_ratio <= 0.0 || temperature_ratio > 1.0 {
        // check valid temperature ratio
        return Err(IsentropicFlowError::InvalidTemperatureRatio);
    }
    let mach_number: f64 = (2.0 * ((1.0 / temperature_ratio) - 1.0)) / (specific_heat_ratio - 1.0).sqrt();
    Ok(mach_number)
}

pub fn calc_mach_from_pressure_ratio(pressure_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if pressure_ratio <= 0.0 || pressure_ratio > 1.0 {
        // check valid pressure ratio
        return Err(IsentropicFlowError::InvalidPressureRatio);
    }
    let mach_number: f64 = (2.0 * ((1.0 / pressure_ratio.powf((specific_heat_ratio - 1.0) / specific_heat_ratio)) - 1.0)) / (specific_heat_ratio - 1.0);
    Ok(mach_number)
}

pub fn calc_mach_from_density_ratio(density_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if density_ratio <= 0.0 || density_ratio > 1.0 {
        // check valid density ratio
        return Err(IsentropicFlowError::InvalidDensityRatio);
    }
    let mach_number: f64 = ((2.0 * (density_ratio.powf((specific_heat_ratio - 1.0) / specific_heat_ratio) - 1.0)) / (specific_heat_ratio - 1.0)).sqrt();
    Ok(mach_number)
}

pub fn calc_mach_from_prandtl_meyer_angle(prandtl_meyer_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    // eta is the value of (m^2 - 1).sqrt()
    let specific_heat_ratio_constant = ((specific_heat_ratio + 1.0) / (specific_heat_ratio - 1.0)).sqrt();
    let prandtl_meyer_function = |eta: f64| {
        specific_heat_ratio_constant * (eta / specific_heat_ratio_constant).atan()
        - eta.atan() - prandtl_meyer_angle    
    };
    let prandtl_meyer_function_1st_derivative = |eta: f64| {
        1.0 / ((eta / specific_heat_ratio_constant).powi(2) + 1.0)
        - 1.0 / (eta.powi(2) + 1.0)
    };
    let eta: f64 = numerics::newton_raphson(
        &prandtl_meyer_function,
        &prandtl_meyer_function_1st_derivative, 
        1.5,
        None,
        None,
    );
    let mach_number: f64 = (eta.powi(2) + 1.0).sqrt();
    Ok(mach_number)
}

pub fn valid_specific_heat_ratio(specific_heat_ratio: f64) -> bool {
    specific_heat_ratio > 1.0
}