#![cfg(test)]

use crate::isentropic::*;
use std::f64::consts::PI;


#[test]
fn test_calc_mach_from_mach_angle() {
    let mach_angle = PI / 6.0;
    let result = calc_mach_from_mach_angle(mach_angle);
    
    match result {
        Ok(mach_number) => {
            println!("mach number is {} for mach angle {}", mach_number, mach_angle);
        },
        Err(error) => {
            panic!("function returned an error: {:?}", error);
        }
    }
}

#[test]
fn test_calc_mach_from_temperature_ratio() {
    let temperature_ratio = 0.8;
    let specific_heat_ratio = 1.4;
    let result = calc_mach_from_temperature_ratio(temperature_ratio, specific_heat_ratio);
    
    match result {
        Ok(mach_number) => {
            println!(
                "mach number is {} for temperature ratio {} and specific heat ratio {}",
                mach_number, temperature_ratio, specific_heat_ratio
            );
        }
        Err(error) => {
            panic!("function returned an error: {:?}", error);
        }
    }
}

#[test]
fn test_calc_mach_from_pressure_ratio() {
    let pressure_ratio = 0.5;
    let specific_heat_ratio = 1.4;
    let result = calc_mach_from_pressure_ratio(pressure_ratio, specific_heat_ratio);
    
    match result {
        Ok(mach_number) => {
            println!(
                "mach number is {} for pressure ratio {} and specific heat ratio {}",
                mach_number, pressure_ratio, specific_heat_ratio
            );
        }
        Err(error) => {
            panic!("function returned an error: {:?}", error);
        }
    }
}

#[test]
fn test_calc_mach_from_density_ratio() {
    let density_ratio = 0.6;
    let specific_heat_ratio = 1.4;
    let result = calc_mach_from_density_ratio(density_ratio, specific_heat_ratio);
    
    match result {
        Ok(mach_number) => {
            dbg!(mach_number);
            println!(
                "mach number is {} for density ratio {} and specific heat ratio {}",
                mach_number, density_ratio, specific_heat_ratio
            );
        }
        Err(error) => {
            panic!("function returned an error: {:?}", error);
        }
    }
}

#[test]
fn test_calc_mach_from_prandtl_meyer_angle() {
    let prandtl_meyer_angle = PI / 4.0;
    let specific_heat_ratio = 1.4;
    let result = calc_mach_from_prandtl_meyer_angle(prandtl_meyer_angle, specific_heat_ratio);

    match result {
        Ok(mach_number) => {
            println!("mach number is: {} for p-m angle of {}", mach_number, prandtl_meyer_angle);
        }
        Err(error) => {
            panic!("function returned an error: {:?}", error);
        }
    }
}