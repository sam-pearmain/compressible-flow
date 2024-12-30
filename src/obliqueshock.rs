#![allow(dead_code)]

use std::result;

use crate::isentropic::{valid_specific_heat_ratio, IsentropicFlowError}; 
use crate::numerics::{bisection, newton_raphson};


pub enum Input {
    UpstreamMach(f64), // value must be given
    DeflectionAngle(f64), // potentially need to change to allow for strong or weak oblique shocks
    ShockAngle(f64),
}

pub enum Output {
    DownstreamMach,
    DeflectionAngle,
    ShockAngle,
    PressureRatio,
    DensityRatio,
    TemperatureRatio,
    StagnationPressureRatio,
    NormalUpstreamMach,
    NormalDownstreamMach,
}

pub struct ObliqueShock {
    upstream_mach: f64,             // M1
    downstream_mach: f64,           // M2
    deflection_angle: f64,          // θ
    shock_angle: f64,               // β
    pressure_ratio: f64,            // p2 / p1
    density_ratio: f64,             // ρ2 / ρ1
    temperature_ratio: f64,         // T2 / T1
    stagnation_pressure_ratio: f64, // p02 / p01
    normal_upstream_mach: f64,      // M1n
    normal_downstream_mach: f64,    // M2n
}

impl ObliqueShock {
    pub fn from_mach_and_shock_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, IsentropicFlowError> {
        // pass
    }
    
    pub fn from_mach_and_deflection_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, IsentropicFlowError> {
        // pass
    }

    pub fn from_mach_and_normal_mach(upstream_mach: f64, normal_upstream_mach: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, IsentropicFlowError> {
        // pass
    }
}

pub fn calculate(input: Vec<Input>, output: Output, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    // pass
}

pub fn calc_downstream_mach() -> Result<f64, IsentropicFlowError> {
    // pass
}

pub fn calc_shock_angle() -> Result<f64, IsentropicFlowError> {
    // pass
}

pub fn calc_deflection_angle() -> Result<f64, IsentropicFlowError> {
    // pass
}

pub fn calc_pressure_ratio() -> Result<f64, IsentropicFlowError> {
    // pass
}

pub fn calc_density_ratio() -> Result<f64, IsentropicFlowError> {
    // pass
}

pub fn calc_temperature_ratio() -> Result<f64, IsentropicFlowError> {
    // pass
}

pub fn calc_max_shock_angle() -> Result<f64, IsentropicFlowError> {
    // calculates the maximum shock angle for a given mach number and specific heat ratio
    // before the oblique shock detatches and becomes a normal bow shock
}

pub fn calc_normal_upstream_mach(upstream_mach: f64, shock_angle: f64) -> Result<f64, IsentropicFlowError> {
    if upstream_mach < 1.0 {
        return Err(IsentropicFlowError::InvalidMachNumber);
    }
    Ok(upstream_mach * shock_angle.sin())
}

pub fn calc_normal_downstream_mach(downstream_mach: f64, shock_angle: f64, deflection_angle: f64) -> Result<f64, IsentropicFlowError> {
    if downstream_mach < 1.0 {
        return Err(IsentropicFlowError::InvalidMachNumber);
    }
    Ok(downstream_mach * (shock_angle - deflection_angle).sin())
}