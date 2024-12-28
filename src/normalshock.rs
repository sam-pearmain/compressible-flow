use crate::numerics;
use crate::isentropic::{ IsentropicFlowError, valid_specific_heat_ratio};


pub enum Input {
    UpstreamMachNumber(f64),
    DownstreamMachNumber(f64),
    TemperatureRatio(f64),
    PressureRatio(f64),
    DensityRatio(f64),
    StagnationTemperatureRatio(f64),
    StagnationPressureRatio(f64),
}

pub enum Output {
    UpstreamMachNumber,
    DownstreamMachNumber,
    TemperatureRatio,
    PressureRatio,
    DensityRatio,
    StagnationTemperatureRatio,
    StagnationPressureRatio,
}

struct NormalShock {
    upstream_mach_number: f64,          // M1
    downstream_mach_number: f64,        // M2
    temperature_ratio: f64,             // T2 / T1
    pressure_ratio: f64,                // p2 / p1
    density_ratio: f64,                 // ρ2 / ρ1
    stagnation_temperature_ratio: f64,  // T02 / T01
    stagnation_pressure_ratio: f64,     // p02 / p01
}

impl NormalShock {
    pub fn new(input: Input, output: Output, specific_heat_ratio: Option<f64>) -> Result<NormalShock, IsentropicFlowError> {
        let specific_heat_ratio = specific_heat_ratio.unwrap_or(1.4);
        if !valid_specific_heat_ratio(specific_heat_ratio) {
            return Err(IsentropicFlowError::InvalidSpecificHeatRatio);
        }
    }
}

pub fn calculate() {

}
