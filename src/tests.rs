#![cfg(test)]

use std::f64::consts::PI;
use crate::isentropic::IsentropicFlow;
use crate::normalshock::NormalShock;
use crate::obliqueshock::ObliqueShock;

#[test]
fn test_isentropic() {
    let mach_number = 2.0;
    let specific_heat_ratio = 1.4;

    match IsentropicFlow::from_mach(mach_number, specific_heat_ratio) {
        Ok(isentropic_flow) => {
            println!("{:?}", isentropic_flow);
        }
        Err(e) => {
            panic!("failed with error: {:?}", e);
        }
    }
}

#[test]
fn test_normal_shock() {
    let upstream_mach = 3.0;
    let specific_heat_ratio = 1.4;

    match NormalShock::from_upstream_mach(upstream_mach, specific_heat_ratio) {
        Ok(normal_shock) => {
            println!("{:?}", normal_shock);
        }
        Err(e) => {
            panic!("failed with error: {:?}", e);
        }
    }
}

#[test]
fn test_oblique_shock() {
    let upstream_mach = 3.0;
    let deflection_angle = PI / 6.0; // 30 deg
    let specific_heat_ratio = 1.4;
    
    match ObliqueShock::from_mach_and_deflection_angle(upstream_mach, deflection_angle, specific_heat_ratio) {
        Ok(oblique_shock) => {
            println!("{:?}", oblique_shock);
        }
        Err(e) => {
            panic!("what the fuckkkkkk {:?}", e);
        }
    }
}

#[test]
fn test_taylor_maccoll() {

}