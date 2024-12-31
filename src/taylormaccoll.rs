#![allow(dead_code)]

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
    shock_turn_angle: f64,                  // ^^^
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