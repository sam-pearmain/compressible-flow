#![allow(dead_code)]

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
    pub fn from_mach_and_deflection_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, IsentropicFlowError> {
        // pass
    }

    pub fn from_mach_and_shock_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, IsentropicFlowError> {
        // pass
    }

    pub fn from_mach_and_normal_mach(upstream_mach: f64, normal_upstream_mach: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, IsentropicFlowError> {
        // pass
    }
}

pub fn calculate(input: Vec<Input>, output: Output, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    // pass
}
