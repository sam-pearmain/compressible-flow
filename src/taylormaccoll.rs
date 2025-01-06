#![allow(dead_code)]

use std::f64::consts::PI;
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
    pub fn from_mach_and_cone_angle(upstream_mach: f64, cone_angle: f64, specific_heat_ratio: f64) -> Result<SupersonicCone, IsentropicFlowError> {
        if cone_angle > PI / 2.0 || cone_angle < 0.0 {
            return Err(IsentropicFlowError::WhatTheFuck);
        }
        // guess a shock angle, use mach angle since mach_angle < shock_angle < PI / 2 rads
        let mut shock_angle = calc_mach_angle_from_mach(upstream_mach)?;

        // need function to calculate cone angle for a given freestream mach and shock angle
        Ok(SupersonicCone)
    }

    pub fn from_mach_and_shock_angle(upstream_mach: f64, shock_angle: f64, specific_heat_ratio: f64) -> Result<SupersonicCone, IsentropicFlowError> {
        let mach_angle: f64 = calc_mach_angle_from_mach(upstream_mach)?;
        if shock_angle > PI / 2.0 || shock_angle - mach_angle < 0.0 {
            return Err(IsentropicFlowError::WhatTheFuck);
        }

        // need function to calculate flow deflection angle 
        Ok(SupersonicCone)
    }

    pub fn from_mach_and_surface_mach() {

    }
}

pub fn solve_taylor_maccoll() {
    
}

pub fn taylor_maccoll(y: (f64, f64), theta: f64, specific_heat_ratio: f64) -> Result<(f64, f64), IsentropicFlowError> {
    if !valid_specific_heat_ratio(specific_heat_ratio) {
        return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
    }
    
    // extract the radial and tangential velocity components
    let radial_velocity: f64 = y.0;
    let tangential_velocity: f64 = y.1;

    // calculate the radial velocity derivative
    let radial_velocity_derivative = tangential_velocity;

    // caclulate the tangential velocity derivative
    let denominator = 
        ((specific_heat_ratio - 1.0) / 2.0) 
        * (1.0 - radial_velocity.powi(2) - tangential_velocity.powi(2)) 
        - tangential_velocity.powi(2);
    let numerator = 
        tangential_velocity.powi(2) * radial_velocity
        - ((specific_heat_ratio - 1.0) / 2.0)
        * (1.0 - radial_velocity.powi(2) - tangential_velocity.powi(2))
        * (2.0 * radial_velocity + tangential_velocity / theta.tan());

    let tangential_velocity_derivative = numerator / denominator;

    Ok((radial_velocity_derivative, tangential_velocity_derivative))
}