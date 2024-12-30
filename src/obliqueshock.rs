#![allow(dead_code)]

use std::f64::consts::PI;
use crate::isentropic::{valid_specific_heat_ratio, IsentropicFlowError}; 
use crate::numerics::bisection;


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
        let deflection_angle: f64 = calc_deflection_angle(upstream_mach, shock_angle, specific_heat_ratio)?;
        let downstream_mach: f64 = calc_downstream_mach(upstream_mach, shock_angle, deflection_angle, specific_heat_ratio)?;
        let pressure_ratio: f64 = calc_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
        let density_ratio: f64 = calc_density_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
        let temperature_ratio: f64 = calc_temperature_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
        let stagnation_pressure_ratio: f64 = pressure_ratio / density_ratio.powf(specific_heat_ratio);
        let normal_upstream_mach: f64 = calc_normal_upstream_mach(upstream_mach, shock_angle)?;
        let normal_downstream_mach: f64 = calc_normal_downstream_mach(downstream_mach, shock_angle, deflection_angle)?;

        Ok(ObliqueShock {
            upstream_mach,
            downstream_mach,
            deflection_angle,
            shock_angle,
            pressure_ratio,
            density_ratio,
            temperature_ratio,
            stagnation_pressure_ratio,
            normal_upstream_mach,
            normal_downstream_mach,
        })
    }
    
    pub fn from_mach_and_deflection_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, IsentropicFlowError> {
        let shock_angle = calc_shock_angle(upstream_mach, deflection_angle, specific_heat_ratio)?;
        ObliqueShock::from_mach_and_shock_angle(upstream_mach, shock_angle, specific_heat_ratio)
    }

    pub fn from_mach_and_normal_mach(upstream_mach: f64, normal_upstream_mach: f64, specific_heat_ratio: f64) -> Result<ObliqueShock, IsentropicFlowError> {
        if upstream_mach <= 1.0 {
            return Err(IsentropicFlowError::InvalidMachNumber);
        }
        if normal_upstream_mach <= 1.0 || normal_upstream_mach > upstream_mach {
            return Err(IsentropicFlowError::InvalidMachNumber);
        }

        let shock_angle = (normal_upstream_mach / upstream_mach).asin();
        ObliqueShock::from_mach_and_shock_angle(upstream_mach, shock_angle, specific_heat_ratio)
    }
}

pub fn calc_downstream_mach(upstream_mach: f64, shock_angle: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    let normal_upstream_mach: f64 = calc_normal_upstream_mach(upstream_mach, shock_angle)?;

    let normal_downstream_mach: f64 = (
        (1.0 + (specific_heat_ratio - 1.0) * normal_upstream_mach.powi(2) / 2.0) /
        (specific_heat_ratio * normal_upstream_mach.powi(2) - (specific_heat_ratio - 1.0) / 2.0)
    ).sqrt();

    let downstream_mach: f64 = normal_downstream_mach / (shock_angle - deflection_angle).cos();
    Ok(downstream_mach)
}

pub fn calc_shock_angle(upstream_mach: f64, deflection_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if upstream_mach <= 1.0 {
        return Err(IsentropicFlowError::InvalidMachNumber);
    }
    let f = |shock_angle: f64| {
        let calculated_deflection_angle = match calc_deflection_angle(upstream_mach, shock_angle, specific_heat_ratio) {
            Ok(value) => value,
            Err(_) => panic!("erm what"),
        };
        return calculated_deflection_angle - deflection_angle
    };

    let lower_bound = deflection_angle;
    let upper_bound = PI / 2.0;

    let shock_angle = bisection(&f, lower_bound, upper_bound, None, None);

    if shock_angle.is_nan() {
        return Err(IsentropicFlowError::MathError);
    }

    Ok(shock_angle)
}

pub fn calc_deflection_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    let tan_deflection_angle = 
        2.0 / shock_angle.tan() * 
        (upstream_mach.powi(2) * shock_angle.sin().powi(2) - 1.0) / 
        (upstream_mach.powi(2) * (specific_heat_ratio + (2.0 * shock_angle).cos()) + 2.0);
    Ok(tan_deflection_angle.atan())
}

pub fn calc_pressure_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    let pressure_ratio: f64 = 
        2.0 * upstream_mach.powi(2) * shock_angle.sin().powi(2) /
        (specific_heat_ratio + 1.0);
    Ok(pressure_ratio)
}

pub fn calc_density_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    let density_ratio: f64 = 
        (specific_heat_ratio + 1.0) * upstream_mach.powi(2) * shock_angle.sin().powi(2) /
        ((specific_heat_ratio - 1.0) * upstream_mach.powi(2) * shock_angle.sin().powi(2) + 2.0);
    Ok(density_ratio)

}

pub fn calc_temperature_ratio(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    let pressure_ratio: f64 = calc_pressure_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
    let density_ratio: f64 = calc_density_ratio(upstream_mach, shock_angle, specific_heat_ratio)?;
    let temperature_ratio: f64 = pressure_ratio * density_ratio;
    Ok(temperature_ratio)
}

pub fn calc_max_shock_angle(upstream_mach: f64, specific_heat_ratio: f64) -> Result<f64, IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    if upstream_mach <= 1.0 {
        return Err(IsentropicFlowError::InvalidMachNumber);
    }

    let sin_max_shock_angle = 
        ((1.0 / (specific_heat_ratio * upstream_mach.powi(2))) * 
        (1.0 +
            ((specific_heat_ratio + 1.0) * (
                (specific_heat_ratio + 1.0) * upstream_mach.powi(4) / 16.0 +
                (specific_heat_ratio - 1.0) * upstream_mach.powi(2) / 2.0 +
                1.0
            ).sqrt())
        )).sqrt();

    if sin_max_shock_angle > 1.0 || sin_max_shock_angle < 0.0 {
        return Err(IsentropicFlowError::MathError);
    }

    let shock_angle: f64 = sin_max_shock_angle.asin();
    Ok(shock_angle)
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