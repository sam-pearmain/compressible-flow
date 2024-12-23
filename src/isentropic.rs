use std::collections::HashMap;
use std::f64::consts::PI;
use roots::{find_root_newton_raphson, SimpleConvergency};

pub enum IsentropicFlowError {
    InvalidMachNumber,
    InvalidMachAngle,
    InvalidTemperatureRatio,
    InvalidPressureRatio,
    InvalidDensityRatio,
    InvalidPrandtlMeyerAngle,
}

pub fn calc_mach_from_mach_angle(mach_angle: f64) -> Result<f64, IsentropicFlowError> {
    if mach_angle < 0.0 || mach_angle > PI / 2.0 {
        // check valid mach angle in radians
        return Err(IsentropicFlowError::InvalidMachAngle)
    }
    let mach_number = 1.0 / mach_angle.sin();
    Ok(mach_number);
}

pub fn calc_mach_from_temperature_ratio(temperature_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {

}

pub fn calc_mach_from_pressure_ratio(pressure_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    
}

pub fn calc_mach_from_density_ratio(density_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    
}

pub fn calc_mach_from_prandtl_meyer_angle(prandtl_meyer_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    
}