use std::time::Instant;
use std::io::Write;
/// This programm generates syntetic, non-equidistant data
/// and produces the Lomb-Scargle periodogram. The behaviour
/// is linear in ofac and hifac and quadratic in N O(ofac * hifac * N^2)


const N: usize = 100;  // number of data points
const T_MAX: f64 = 100.;  // recorded seconds
const FREQUENCY: f64 = 2.*PI*0.55;  // frequency of the syntetic signal

const PI: f64 = std::f64::consts::PI;


#[allow(dead_code)]
fn read_from_file(_filename: &str, _abs: i32, _ord: i32) -> (Vec<f64>, Vec<f64>) {
  let mut _abscissa: Vec<f64> = Vec::new();
  let mut _ordinate: Vec<f64> = Vec::new();

  (_abscissa, _ordinate)
}


fn write_timedata_to_file(time: [f64; N], amplitude: [f64; N], filename: &str) -> std::io::Result<()> {
  if time.len() != amplitude.len() {
    panic!("The arrays, containing the temporal data, to write into a datafile have unequal lengths.")
  }
  let mut file = std::fs::OpenOptions::new()
                                      .create(true)
                                      .write(true)
                                      .truncate(true)
                                      .open(filename)?;
  file.write_all(b"# time [s]; amplitude [arb. u]\n")
        .expect("Couldn't write to file.");
  for i in 0..N {
    let buffer = format!("{:.4}; {:.4}\n", time[i], amplitude[i]);
    file.write_all(buffer.as_bytes())
        .expect("Couldn't write to file.");
    file.flush()?;
  }
  Ok(())
}


fn write_specdata_to_file(freq: &Vec<f64>, spec: &Vec<f64>, filename: &str) -> std::io::Result<()> {
  if freq.len() != spec.len() {
    panic!("The arrays, containing the spectal data, to write into a datafile have unequal lengths.")
  }
  let mut file = std::fs::OpenOptions::new()
                                      .create(true)
                                      .write(true)
                                      .truncate(true)
                                      .open(filename)?;
  file.write_all(b"# fequency [1/s]; amplitude [arb. u]\n")
        .expect("Couldn't write to file.");
  for i in 0..freq.len() {
    let buffer = format!("{:.4}; {:.4}\n", freq[i], spec[i]);
    file.write_all(buffer.as_bytes())
        .expect("Couldn't write to file.");
    file.flush()?;
  }
  Ok(())
}


fn random_times() -> [f64; N] {
  // create an array of random numbers
  let mut buffer: [f64; N] = [0.; N];
  for item in buffer.iter_mut().take(N) {
      *item = T_MAX * rand::random::<f64>();
    }

  // sort that so the values are monotonically increasing
  buffer.sort_by(|a, b| a.partial_cmp(b).unwrap());

  buffer
}


fn syntetic_signal(mut payload: [f64; N]) -> [f64; N] {
  for item in payload.iter_mut().take(N) {
      *item *= FREQUENCY;
      *item = 1.0 * f64::sin(*item);  //+ 0.66 * f64::sin(0.8 * *item);
  }
  payload
}


fn mean_and_variance(payload: &[f64; N]) -> (f64, f64) {
  let mut mean: f64 = 0.;
  for x in payload {
    mean += x;
  }
  mean /= payload.len() as f64;

  let mut var: f64 = 0.;
  for x in payload {
    var += f64::powi(x - mean, 2);
  }
  var /= payload.len() as f64;

  (mean, var)
}


fn return_max_and_min(payload: &[f64; N]) -> (f64, f64) {
  let mut max: f64 = payload[0];
  let mut min: f64 = payload[0];
  for item in payload.iter().take(N) {
    if *item > max {max = *item;}
    if *item < min {min = *item;}
  }
  (max, min)
}


#[allow(clippy::too_many_arguments)]
fn lomb_spectrum(time: &[f64; N], amplitude: &[f64; N], ofac: &f64, hifac: &f64,
                 freq: &mut Vec<f64>, spec: &mut Vec<f64>,
                 nout: &mut usize, jmax: &mut usize, prob: &mut f64) {
  const TWO_PI: f64 = std::f64::consts::TAU;  // equal to 2*PI
  *nout = (0.5*ofac*hifac*N as f64) as usize;

  // resize freq and spec to nout
  freq.resize(*nout, 0.);
  spec.resize(*nout, 0.);

  // get mean and variance of the input data (amplitude)
  let (mean, var): (f64, f64) = mean_and_variance(amplitude);

  // throw an error if var == 0, for a variance of 0 leads to a bad state.
  if var == 0. {
    panic!("Variance is 0.")
  }

  // go throug the data to get the range of abscissas
  let (tmax, tmin): (f64, f64) = return_max_and_min(time);

  // calculate values for the trigonometric recurrences
  let time_diff: f64 = tmax - tmin;
  let time_avr: f64 = 0.5 * (tmax + tmin);
  let mut spec_max: f64 = 0.;
  let mut freq_now: f64 = 1.0 / (time_diff * ofac);
  let mut wpr: [f64; N] = [0.; N];
  let mut wpi: [f64; N] = [0.; N];
  let mut wr:  [f64; N] = [0.; N];
  let mut wi:  [f64; N] = [0.; N];
  for i in 0..N {
    let arg: f64 = TWO_PI * ((time[i] - time_avr) * freq_now);
    //wpr[i] = -2.0 * f64::sqrt(f64::sin(0.5*arg));
    wpr[i] = -2.0 * f64::powi(f64::sin(0.5*arg), 2);
    wpi[i] = f64::sin(arg);
    wr[i]  = f64::cos(arg);
    wi[i]  = wpi[i];
  }

  // main loop over the frequencys
  for i in 0..*nout {
    freq[i] = freq_now;
    let mut sumsh: f64 = 0.;
    let mut sumc:  f64 = 0.;

    // loop over the data to get tau
    for i in 0..N {
      let c: f64 = wr[i];
      let s: f64 = wi[i];
      sumsh += s*c;
      sumc  += (c-s)*(c+s);
    }
    let wtau: f64 = 0.5 * f64::atan2(2.*sumsh, sumc);
    let swtau: f64 = f64::sin(wtau);
    let cwtau: f64 = f64::cos(wtau);

    // loop over the data to get the spec value
    let mut sums: f64 = 0.;
    let mut sumc: f64 = 0.;
    let mut sumsy: f64 = 0.;
    let mut sumcy: f64 = 0.;
    for i in 0..N {
      let c: f64 = wr[i];
      let s: f64 = wi[i];
      let ss: f64 = s * cwtau - c * swtau;
      let cc: f64 = c * cwtau + s * swtau;
      sums += ss * ss;
      sumc += cc * cc;
      let yy: f64 = amplitude[i] - mean;
      sumsy += yy * ss;
      sumcy += yy * cc;
      let wtemp:f64 = wr[i];
      wr[i] += wr[i] * wpr[i] - wi[i] * wpi[i];
      wi[i] += wi[i] * wpr[i] + wtemp * wpi[i];
    }
    spec[i] = 0.5 * (sumcy * sumcy / sumc + sumsy * sumsy / sums) / var;
    if spec[i] > spec_max {
      spec_max = spec[i];
      *jmax = i;
    }

    freq_now += 1.0 / (ofac * time_diff);
  }

  // evaluate the statistical significance of the maximum
  let expy: f64 = f64::exp(-spec_max);
  let eff_m: f64 = 2.0 * *nout as f64 / ofac;
  *prob = eff_m * expy;
  if *prob > 0.01 {
    *prob = 1.0 - f64::powf(1.0 - expy, eff_m);
  }
}


fn give_significance_levels(level: f64, nout: &usize, ofac: &f64) -> f64 {
  let eff_m: f64 = 2.0 * *nout as f64 / ofac;

  -f64::ln(1. - f64::powf(1. - level, 1. / eff_m))
}


fn main() -> Result<(), std::io::Error> {
  let now = Instant::now();

  // generate a random time signal
  let time: [f64; N] = random_times();

  // generate a syntetic amplitude signal
  let amplitude: [f64; N] = syntetic_signal(time);

  let ofac: f64 = 4.;  // oversamplig factor (normally: 4)
  let hifac: f64 = 2.;  // factor for the highest frequency: freq_max / freq_nyquist
  let mut freq: Vec<f64> = Vec::new();  // output vector for the frequencys
  let mut spec: Vec<f64> = Vec::new();  // output vector for the spectrum
  let mut freq_len: usize = 0;
  let mut jmax: usize = 0;
  let mut prob:f64 = 0.;
  lomb_spectrum(&time, &amplitude, &ofac, &hifac, &mut freq,
                &mut spec, &mut freq_len, &mut jmax, &mut prob);

  write_timedata_to_file(time, amplitude, "data.temp")?;
  write_specdata_to_file(&freq, &spec, "data.spec")?;

  //println!("freq_len: {}", freq_len);
  //println!("jmax: {}", jmax);
  println!("Probability of the highest peak to be random noise: {:.0e}", prob);

  println!("Significance level of {{}} at {{}}");
  println!("0.0001 at {:.2}", give_significance_levels(0.0001, &freq_len, &ofac));
  println!(" 0.001 at {:.2}", give_significance_levels( 0.001, &freq_len, &ofac));
  println!("  0.01 at {:.2}", give_significance_levels(  0.01, &freq_len, &ofac));
  println!("   0.1 at {:.2}", give_significance_levels(   0.1, &freq_len, &ofac));
  
  let elapsed = now.elapsed();
  println!("{:.2?}", elapsed);

  Ok(())
}
