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
    upstream_mach: f64,
    downstream_mach: f64, 
    deflection_angle: f64,
    shock_angle: f64,
    pressure_ratio: f64,
    density_ratio: f64,
    temperature_ratio: f64,
    stagnation_pressure_ratio: f64,
    normal_upstream_mach: f64,
    normal_downstream_mach: f64,
}

