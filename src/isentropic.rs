use std::f64::consts::PI;
use crate::numerics::*;


#[derive(Debug)]
pub enum IsentropicFlowError {
    InvalidMachNumber,
    InvalidMachAngle,
    InvalidTemperatureRatio,
    InvalidPressureRatio,
    InvalidDensityRatio,
    InvalidPrandtlMeyerAngle,
    InvalidSpecificHeatRatio,
    WhatTheFuck,
    MathError,
}

pub enum Input {
    MachNumber(f64),
    MachAngle(f64),
    TemperatureRatio(f64),
    PressureRatio(f64),
    DensityRatio(f64),
    PrandtlMeyerAngle(f64),
}

pub enum Output {
    MachNumber,
    MachAngle,
    TemperatureRatio,
    PressureRatio,
    DensityRatio,
    PrandtlMeyerAngle,
}

#[derive(Debug)]
pub struct IsentropicFlow {
    mach_number: f64,           // M
    mach_angle: f64,            // μ
    temperature_ratio: f64,     // T / T0
    pressure_ratio: f64,        // p / p0
    density_ratio: f64,         // ρ / ρ0
    prandtl_meyer_angle: f64,   // 𝒱(M)
}

impl IsentropicFlow {
    pub fn from_mach(mach_number: f64, specific_heat_ratio: f64) -> Result<IsentropicFlow, IsentropicFlowError> {
        if !valid_specific_heat_ratio(specific_heat_ratio) {
            return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
        }
        if mach_number < 0.0 {
            return Err(IsentropicFlowError::InvalidMachNumber);
        }
        
        let mach_angle = calc_mach_angle_from_mach(mach_number)?;
        let temperature_ratio = calc_temperature_ratio_from_mach(mach_number, specific_heat_ratio)?;
        let pressure_ratio = calc_pressure_ratio_from_mach(mach_number, specific_heat_ratio)?;
        let density_ratio = calc_density_ratio_from_mach(mach_number, specific_heat_ratio)?;
        let prandtl_meyer_angle = prandtl_meyer_function(mach_number, specific_heat_ratio)?;

        Ok(IsentropicFlow{
            mach_number,
            mach_angle,
            temperature_ratio,
            pressure_ratio,
            density_ratio,
            prandtl_meyer_angle,
        })
    }
}

pub fn calculate(output: Output, input: Input, specific_heat_ratio: Option<f64>) -> Result<f64, IsentropicFlowError> {
    // simple calculator function for if you're lazy 
    let specific_heat_ratio = specific_heat_ratio.unwrap_or(1.4);
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }

    let isentropic = match input {
        Input::MachNumber(value) => {
            IsentropicFlow::from_mach(value, specific_heat_ratio)?
        }
        Input::MachAngle(value) => {
            let mach_number = calc_mach_from_mach_angle(value)?;
            IsentropicFlow::from_mach(mach_number, specific_heat_ratio)?
        }
        Input::TemperatureRatio(value) => {
            let mach_number = calc_mach_from_temperature_ratio(value, specific_heat_ratio)?;
            IsentropicFlow::from_mach(mach_number, specific_heat_ratio)?
        }
        Input::PressureRatio(value) => {
            let mach_number = calc_mach_from_pressure_ratio(value, specific_heat_ratio)?;
            IsentropicFlow::from_mach(mach_number, specific_heat_ratio)?
        }
        Input::DensityRatio(value) => {
            let mach_number = calc_mach_from_density_ratio(value, specific_heat_ratio)?;
            IsentropicFlow::from_mach(mach_number, specific_heat_ratio)?
        }
        Input::PrandtlMeyerAngle(value) => {
            let mach_number = calc_mach_from_prandtl_meyer_angle(value, specific_heat_ratio)?;
            IsentropicFlow::from_mach(mach_number, specific_heat_ratio)?
        }
    };

    match output {
        Output::MachNumber => {Ok(isentropic.mach_number)}
        Output::MachAngle => {Ok(isentropic.mach_angle)}
        Output::TemperatureRatio => {Ok(isentropic.temperature_ratio)}
        Output::PressureRatio => {Ok(isentropic.pressure_ratio)}
        Output::DensityRatio => {Ok(isentropic.density_ratio)}
        Output::PrandtlMeyerAngle => {Ok(isentropic.prandtl_meyer_angle)}
    }
}

pub fn calc_mach_angle_from_mach(mach_number: f64) -> Result<f64, IsentropicFlowError> {
    if mach_number < 0.0 {
        return Err(IsentropicFlowError::InvalidMachNumber);
    }
    let mach_angle = (1.0 / mach_number).asin();
    Ok(mach_angle)
}

pub fn calc_pressure_ratio_from_mach(mach_number: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let pressure_ratio: f64 = (1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2)).powf(-specific_heat_ratio / (specific_heat_ratio - 1.0));
    Ok(pressure_ratio)
}

pub fn calc_temperature_ratio_from_mach(mach_number: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let temperature_ratio: f64 = (1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2)).powi(-1);
    Ok(temperature_ratio)
}

pub fn calc_density_ratio_from_mach(mach_number: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let density_ratio: f64 = (1.0 + (specific_heat_ratio - 1.0) / 2.0 * mach_number.powi(2)).powf(-1.0 / (specific_heat_ratio - 1.0));
    Ok(density_ratio)
}

pub fn calc_mach_from_speed_of_sound(velocity: f64, speed_of_sound: f64) -> Result<f64, IsentropicFlowError> {
    Ok(velocity / speed_of_sound)
}

pub fn prandtl_meyer_function(mach_number: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    if mach_number <= 1.0 {
        return Err(IsentropicFlowError::InvalidMachNumber);
    }
    let gamma_ratio = (specific_heat_ratio - 1.0) / (specific_heat_ratio + 1.0);
    let sqrt_gamma_ratio = gamma_ratio.sqrt();
    let prandtl_meyer_angle: f64 = (
        (1.0 / sqrt_gamma_ratio) // 1st term
        * (sqrt_gamma_ratio * (mach_number.powi(2) - 1.0).sqrt()).atan()) // 2nd term
        - (mach_number.powi(2) - 1.0).sqrt().atan(); // 3rd term
    Ok(prandtl_meyer_angle)
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
    let mach_number: f64 = (2.0 * ((1.0 / temperature_ratio) - 1.0) / (specific_heat_ratio - 1.0)).sqrt();
    Ok(mach_number)
}

pub fn calc_mach_from_pressure_ratio(pressure_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    if pressure_ratio <= 0.0 || pressure_ratio > 1.0 {
        // check valid pressure ratio
        return Err(IsentropicFlowError::InvalidPressureRatio);
    }
    let mach_number: f64 = (2.0 * ((1.0 / pressure_ratio.powf((specific_heat_ratio - 1.0) / specific_heat_ratio)) - 1.0) / (specific_heat_ratio - 1.0)).sqrt();
    Ok(mach_number)
}

pub fn calc_mach_from_density_ratio(density_ratio: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    if density_ratio <= 0.0 || density_ratio > 1.0 {
        // check valid density ratio
        return Err(IsentropicFlowError::InvalidDensityRatio);
    }
    let mach_number: f64 = ((2.0 * ((1.0 / density_ratio.powf(specific_heat_ratio - 1.0)) - 1.0)) / (specific_heat_ratio - 1.0)).sqrt();
    Ok(mach_number)
}

pub fn calc_mach_from_prandtl_meyer_angle(prandtl_meyer_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    // eta is the value of (m^2 - 1).sqrt()
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    let alpha = ((specific_heat_ratio + 1.0) / (specific_heat_ratio - 1.0)).sqrt(); // just a constant to make things easier
    let f = |eta: f64| {
        alpha * (eta / alpha).atan()
        - eta.atan() - prandtl_meyer_angle    
    };
    let df = |eta: f64| {
        1.0 / ((eta / alpha).powi(2) + 1.0)
        - 1.0 / (eta.powi(2) + 1.0)
    };
    let eta: f64 = newton_raphson(&f, &df, 1.5, None, None);
    let mach_number: f64 = (eta.powi(2) + 1.0).sqrt();
    Ok(mach_number)
}

pub fn valid_specific_heat_ratio(specific_heat_ratio: f64) -> bool {
    // specific heat ratio must be greater than 1
    specific_heat_ratio > 1.0
}