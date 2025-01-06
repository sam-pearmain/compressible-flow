#![cfg(test)]

use crate::isentropic::IsentropicFlow;
use crate::normalshock::NormalShock;


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

}

#[test]
fn test_taylor_maccoll() {

}