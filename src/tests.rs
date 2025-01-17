#![cfg(test)]

use std::f64::consts::PI;
use crate::isentropic::IsentropicFlow;
use crate::normalshock::NormalShock;
use crate::obliqueshock::{self, ObliqueShock};
use crate::taylormaccoll;

#[test]
fn test_isentropic() {
    // working
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
    // working
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
fn test_calc_shock_angle_from_deflection_angle() {
    // working
    let upstream_mach = 3.0;
    let deflection_angle = PI / 6.0; // 30 deg
    let specific_heat_ratio = 1.4;
    
    match obliqueshock::calc_shock_angle(upstream_mach, deflection_angle, specific_heat_ratio) {
        Ok(shock_angle) => {println!("shock angle: {:?}", shock_angle);}
        Err(e) => {panic!("what the fuckkkkkk {:?}", e);}
    }
}

#[test]
fn test_calc_downstream_mach_number() {
    // passed
    let upstream_mach = 3.0;
    let shock_angle = PI / 4.0;
    let specific_heat_ratio = 1.4;

    match obliqueshock::calc_downstream_mach_from_shock_angle(upstream_mach, shock_angle, specific_heat_ratio) {
        Ok(downstream_mach) => {println!("downstream mach: {:?}", downstream_mach);}
        Err(e) => {panic!("what now: {:?}", e);}
    }
}

#[test]
fn test_oblique_shock() {
    // passed 
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
fn test_solve_taylor_maccoll() {
    // yeah it doesnt quite work
    // starting stuff
    use std::fs::File;
    use std::io::Write;

    let mach_number = 4.0;
    let shock_angle = PI / 6.0; // 30 deg
    let specific_heat_ratio = 1.4;
    
    // get deflection angle
    let deflection_angle = obliqueshock::calc_deflection_angle(mach_number, shock_angle, specific_heat_ratio).expect("what");

    // get downstream velocity components
    let downstream_mach = obliqueshock::calc_downstream_mach_from_shock_angle(mach_number, shock_angle, specific_heat_ratio).expect("erm");
    let radial_downstream_mach = downstream_mach * (shock_angle - deflection_angle).cos();
    let tangential_downstream_mach = - downstream_mach *(shock_angle - deflection_angle).sin();
    
    assert!((downstream_mach - (tangential_downstream_mach.powi(2) + radial_downstream_mach.powi(2)).sqrt()).abs() < 1e-6);

    let mut file = File::create("results.txt").expect("failed");

    match taylormaccoll::solve_taylor_maccoll(
        (radial_downstream_mach, tangential_downstream_mach), 
        shock_angle, 
        0.0, 
        specific_heat_ratio,
        true, 
        Some(20000),
    ) {
        Ok((velocity_components, thetas)) => {
            assert_eq!(velocity_components.len(), thetas.len());
            for ((radial, tangential), theta) in velocity_components.iter().zip(thetas.iter()) {
                writeln!(
                    file, 
                    "velocity components: ({:.4}, {:.4}), theta: {:.4}",
                    radial, tangential, theta
                ).expect("failed");
            }
        }
        Err(e) => {
            panic!("please no: {:?}", e);
        }
    }
}

#[test]
fn test_taylor_maccoll() {
    // something wrong with my ratios
    let mach_number = 4.0;
    let shock_angle = PI / 6.0; // 30 deg
    let specific_heat_ratio = 1.4;

    match taylormaccoll::SupersonicCone::from_mach_and_shock_angle(mach_number, shock_angle, specific_heat_ratio) {
        Ok(supersonic_cone) => {
            println!("{:?}", supersonic_cone);
        }
        Err(e) => {
            panic!("why, {:?}", e)
        }
    }
}