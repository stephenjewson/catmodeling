use extendr_api::prelude::*;

/// Return string `"Hello world!"` to R.
/// Just a test. 
/// @return A message to the screen.
#[extendr]
fn rust_hello_world() -> &'static str {
    "Hello world!"
}

/// Square an integer.
/// Just an example.
/// @param x An integer.
/// @return The square of the input.
#[extendr]
fn rust_square(x: i32) -> i32 {
    x*x
}

/// Adds two vectors element-wise up to length n.
/// 
/// @param x A reference to the first input vector (slice of f64)
/// @param y A reference to the second input vector (slice of f64)
/// @param n The number of elements to process
/// @return The sum of the two input vectors, as a vector.
#[extendr]
pub fn addvectors(x: &[f64], y: &[f64], n: usize) -> Vec<f64> {
    // z=numeric(n)
    // In R, numeric(n) creates a double-precision vector of length n initialized to 0.
    // In Rust, vec![0.0; n] creates a Vec<f64> of length n initialized to 0.0.
    let mut z1 = vec![0.0; n];
    let mut z2 = vec![0.0; n];

    for i in 0..n {
        z1[i] = x[i] + y[i];
        z2[i] = x[i] - y[i];
    }

    // return(z)
    return z1;
}

/// Function to calculate the Haversine distance between two points.
/// This matches the functionality of geosphere::distHaversine in R.
/// @param lon1 Lon of first point.
/// @param lat1 Lat of first point.
/// @param lon2 Lon of second point.
/// @param lat2 You guess.
/// @return Returns distance in meters.
/// @return A message to the screen.
#[extendr]
fn rust_dist_haversine(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let r = 6378137.0; // Radius of Earth in meters (geosphere default semi-major axis)
    let phi1 = lat1.to_radians();
    let phi2 = lat2.to_radians();
    let lam1 = lon1.to_radians();
    let lam2 = lon2.to_radians();

    let d_phi = phi2 - phi1;
    let d_lam = lam2 - lam1;

    let a = (d_phi / 2.0).sin().powi(2)
        + phi1.cos() * phi2.cos() * (d_lam / 2.0).sin().powi(2);
    
    // Equivalent to 2 * asin(sqrt(a)) but using atan2 for consistency with geosphere
    let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());

    r * c
}

/// Hours clause function
/// For a given start event, calculates all possible hours clauses, returns new losses
/// @param hcdays           Definition of the hours clause in days
/// @param hckm             Definition of the hours clause in km
/// @param temploss         The losses by event for this year  
/// @param thisyearday      The days by event for this year
/// @param thisyearlon      The lons by event for this year   
/// @param thisyearlat      The lats by event for this year
/// @param nevents_in_year  The number of events in the full year
/// @param k                The starting event
/// @return Returns a vector of the new losses
#[extendr]
pub fn rust_hours_clause_apply_part1(
    hcdays: f64,
    hckm: f64,
    mut temploss: Vec<f64>,
    thisyearday: &[f64],
    thisyearlon: &[f64],
    thisyearlat: &[f64],
    nevents_in_year: usize,
    k: usize,
) -> Vec<f64> {

    // apply the hours clause to event k in the year

    // now test whether each subsequent event is captured within the hours clause
    for kk in (k + 1)..=nevents_in_year {
    //for kk in (k)..(nevents_in_year) {

        // hcdays
        let daydiff = (thisyearday[kk-1] - thisyearday[k-1]).abs();

        // space
        let lon1 = thisyearlon[kk-1];
        let lat1 = thisyearlat[kk-1];
        let lon2 = thisyearlon[k-1];
        let lat2 = thisyearlat[k-1];
        // #					cat("lons=",lon1,lon2,"\n")
        // #					cat("lats=",lat1,lat2,"\n")
        
        // kmdiff=0.001*distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
        let kmdiff = 0.001 * rust_dist_haversine(lon1, lat1, lon2, lat2);
        // #					cat("kmdiff=",kmdiff,"\n")

        // if it falls within the hours clause, move the loss onto the earlier event
        if (daydiff <= hcdays) && (kmdiff <= hckm) {
            temploss[k-1] = temploss[k-1] + temploss[kk-1];
            temploss[kk-1] = 0.0;
        }
    }

    return temploss;
    
}

/// Hours clause function
/// For a given start event, calculates all possible hours clauses, returns new losses AND a counter
/// @param hcdays           Definition of the hours clause in days
/// @param hckm             Definition of the hours clause in km
/// @param temploss         The losses by event for this year  
/// @param thisyearday      The days by event for this year
/// @param thisyearlon      The lons by event for this year   
/// @param thisyearlat      The lats by event fo rthis year
/// @param nevents_in_year  The number of events in the full year
/// @param k                The starting event
/// @return Returns a vector of the new losses
#[extendr]
pub fn rust_hours_clause_apply_part2(
    hcdays: f64,
    hckm: f64,
    mut temploss: Vec<f64>,
    thisyearday: &[f64],
    thisyearlon: &[f64],
    thisyearlat: &[f64],
    nevents_in_year: usize,
    k: usize,
) -> extendr_api::Robj {

    // #	cat("inside applyhck\n")
    // #	cat("hcdays=",hcdays,"\n")
    // #	cat("temploss=",temploss,"\n")
    // #	cat("thisyearday=",thisyearday,"\n")
    // #	cat("nevents_in_year=",nevents_in_year,"\n")
    // #	cat("k=",k,"\n")

    // counter for number of events making up a loss
    let mut hccounter = vec![0.0; temploss.len()];

    // apply the hours clause to event k in the year

    // now test whether each subsequent event is captured within the hours clause
//    for kk in (k + 1)..nevents_in_year {
    for kk in (k+1)..=(nevents_in_year) {

        // hcdays
        let daydiff = (thisyearday[kk-1] - thisyearday[k-1]).abs();

        // space
        let lon1 = thisyearlon[kk-1];
        let lat1 = thisyearlat[kk-1];
        let lon2 = thisyearlon[k-1];
        let lat2 = thisyearlat[k-1];
        // #					cat("lons=",lon1,lon2,"\n")
        // #					cat("lats=",lat1,lat2,"\n")
        
        // kmdiff=0.001*distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
        let kmdiff = 0.001 * rust_dist_haversine(lon1, lat1, lon2, lat2);
        // #					cat("kmdiff=",kmdiff,"\n")

        // if it falls within the hours clause, move the loss onto the earlier event
        if (daydiff <= hcdays) && (kmdiff <= hckm) {
            temploss[k-1] = temploss[k-1] + temploss[kk-1];
            temploss[kk-1] = 0.0;
            hccounter[k-1] = hccounter[k-1] + 1.0;
            hccounter[kk-1] = -1.0;
        }
    }

    list!(temploss=temploss,hccounter=hccounter).into()
    
}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
   mod catmodeling;
   fn rust_hello_world;
   fn rust_square;
   fn rust_hours_clause_apply_part1;
   fn rust_hours_clause_apply_part2;
   fn rust_dist_haversine;
   fn addvectors;
}

