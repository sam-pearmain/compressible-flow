pub mod isentropic;
pub mod numerics;


pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use isentropic::calc_mach_from_prandtl_meyer_angle;

    use super::*;

    #[test]
    fn test_calc_mach_from_prandtl_meyer_angle() {
        let prandtl_meyer_angle = std::f64::consts::PI / 4.0;
        let specific_heat_ratio = 1.4;
        let mach_number_result = calc_mach_from_prandtl_meyer_angle(prandtl_meyer_angle, specific_heat_ratio);
    
        match mach_number_result {
            Ok(mach_number) => {
                dbg!(mach_number); // Debug prints the Mach number with additional context
                println!("mach number is: {} for p-m angle of {}", mach_number, prandtl_meyer_angle);
            }
            Err(error) => {
                panic!("function returned an error: {:?}", error);
            }
        }
    }
}
