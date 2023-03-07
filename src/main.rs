use arrayvec::ArrayVec;
use std::io::Write;

const N: usize = 100;  // number of data points
const T_MAX: f64 = 10.;  // recorded seconds

const PI: f64 = std::f64::consts::PI;
const FREQUENCY: f64 = 2.*PI*0.25;

fn write_to_file(time: [f64; N], amplitude: [f64; N]) -> std::io::Result<()> {
  let mut file = std::fs::OpenOptions::new()
                                      .create(true)
                                      .write(true)
                                      .truncate(true)
                                      .open("data.dat")?;
  for i in 0..N {
    let buffer = format!("{:.4}; {:.4}\n", time[i], amplitude[i]);
    file.write_all(buffer.as_bytes())
        .expect("Couldn't write to file.");  // ?
    file.flush()?;
  }
  Ok(())
}

fn signal(mut time: [f64; N]) -> [f64; N] {
  for item in time.iter_mut().take(N) {
      *item *= FREQUENCY;
      *item = f64::sin(*item);
  }
  time
}

fn main() -> Result<(), std::io::Error> {
  // create an array of random numbers
  let mut array_vec: ArrayVec<f64, N> = ArrayVec::new();
  for _ in 0..N {
      array_vec.push(T_MAX * rand::random::<f64>());
    }

  // put that into an array and sort it
  let mut time: [f64; N] = array_vec.into_inner().unwrap();
  time.sort_by(|a, b| a.partial_cmp(b).unwrap());

  // generate a syntetic signal
  let amplitude: [f64; N] = signal(time);

  write_to_file(time, amplitude)?;

  Ok(())
}
