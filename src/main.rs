use std::{
    fs::OpenOptions,
    io::{BufWriter, Write},
};

//czarna magia

struct BoxedFunction {
    f: Box<dyn Fn(&[f64; 15], &f64, &f64, &f64) -> f64>,
}

impl BoxedFunction {
    fn new<F>(f: F) -> BoxedFunction
    where
        F: Fn(&[f64; 15], &f64, &f64, &f64) -> f64 + 'static,
    {
        BoxedFunction { f: Box::new(f) }
    }
}
fn pushinator(
    eq: impl Fn(&[f64; 15], &f64, &f64, &f64) -> f64 + 'static,
    eqs: &mut Vec<BoxedFunction>,
) {
    eqs.push(BoxedFunction::new(eq));
}
macro_rules! zip {
    ($x: expr) => ($x);
    ($x: expr, $($y: expr), +) => (
        $x.iter().zip(
            zip!($($y), +))
    )
}

fn main() {
    let _path1 = "/home/kartonrealista/actual_code/praca_mgr_symulacja/ptau1000.csv";
    let _path2 = "/home/kartonrealista/actual_code/praca_mgr_symulacja/stezenia.csv";
    let path1win =
        r"C:\Users\admin\Desktop\MTHOMAS\x\fenantrolina9zm\praca_mgr_symulacja\ptau4000.csv";
    let path2win =
        r"C:\Users\admin\Desktop\MTHOMAS\x\fenantrolina9zm\praca_mgr_symulacja\stezenia4.csv";

    let mut h = 1e-3;
    let mut t = 0.0;

    //stałe
    let k_edta = 10.0_f64.powf(9.0);
    let k_oh2: f64 = 1.46 * 10.0_f64.powf(-2.0);
    let k_edta2: f64 = 30.0;
    let k_prim_oxid: f64 = 5000.0;
    let k_m5: f64 = 7.5e-4;
    let k_m7: f64 = 500.0;
    let k_m8: f64 = 1000.0;
    let k_m9: f64 = 1.0;
    let k_m18: f64 = 3000.0;
    let k_m19: f64 = 2.0e6;
    let k_m20: f64 = 10.0_f64.powf(5.0);
    let k_6prim: f64 = 6.93 * 10.0_f64.powf(-3.0);
    let k_8prim: f64 = 5.0 * 10.0_f64.powf(4.0);
    let k_m21: f64 = 2000.0;
    let k_m3: f64 = 1e-2;
    let k_m4: f64 = 20.0;
    let k_cu21: f64 = 100.0;
    let k_phen1 = 8e12;
    let k_phen2 = 7e-6;
    //równania
    let eq1 = move |c: &[f64; 15], _c_h2o2: &f64, c_ho2min: &f64, _c_scn: &f64| {
        k_6prim * c[3] - k_8prim * c[1] * c_ho2min
    };
    let eq2 = move |c: &[f64; 15], _c_h2o2: &f64, c_ho2min: &f64, _c_scn: &f64| {
        k_m20 * c[4] * c[9] + k_8prim * c[1] * c_ho2min
            - k_cu21 * c[2] * c_ho2min
            - k_edta * c_ho2min.powf(2.0) * c[2] * c[11]
            + k_prim_oxid * c[9] * (c[12] + c[13])
    };
    let eq3 = move |c: &[f64; 15], _c_h2o2: &f64, c_ho2min: &f64, _c_scn: &f64| {
        k_cu21 * c[2] * c_ho2min - k_6prim * c[3] - k_m3 * c[3]
    };
    let eq4 = move |c: &[f64; 15], _c_h2o2: &f64, _c_ho2min: &f64, _c_scn: &f64| {
        k_m3 * c[3] - k_m20 * c[4] * c[9]
    };
    let eq5 = move |c: &[f64; 15], c_h2o2: &f64, _c_ho2min: &f64, _c_scn: &f64| {
        k_6prim * c[3] - k_m21 * c[5] * c[8] + k_m3 * c[3] - 2.0 * k_m4 * c[5] * c[5]
            + 2.0 * k_oh2 * c_h2o2 * c[12]
            - 2.0 * k_phen1 * c[5] * c[5] * c[14]

    };
    let eq6 = move |c: &[f64; 15], c_h2o2: &f64, _c_ho2min: &f64, c_scn: &f64| {
        k_m5 * c_h2o2 * c_scn - k_m7 * c[6] * c[6]
    };
    let eq7 = move |c: &[f64; 15], _c_h2o2: &f64, _c_ho2min: &f64, _c_scn: &f64| {
        k_m7 * c[6] * c[6] - k_m18 * c[7] * c[8] + k_m19 * c[9] * c[9]
            - 2.0 * k_m9 * c[7] * c[7]
            - k_m8 * c[6] * c[7]
    };
    let eq8 = move |c: &[f64; 15], _c_h2o2: &f64, _c_ho2min: &f64, _c_scn: &f64| {
        k_m9 * c[7] * c[7] - k_m18 * c[7] * c[8] + k_m19 * c[9] * c[9] + k_m20 * c[4] * c[9]
            - k_m21 * c[5] * c[8]
            + k_prim_oxid * c[9] * (c[12] + c[13])
    };
    let eq9 = move |c: &[f64; 15], _c_h2o2: &f64, _c_ho2min: &f64, _c_scn: &f64| {
        2.0 * k_m18 * c[7] * c[8]
            - 2.0 * k_m19 * c[9] * c[9]
            - k_m20 * c[4] * c[9]
            - k_prim_oxid * c[9] * (c[12] + c[13])
            
    };
    let eq10 = move |c: &[f64; 15], _c_h2o2: &f64, _c_ho2min: &f64, _c_scn: &f64| {
        k_m4 * c[5] * c[5] + k_prim_oxid * c[2] * c[9]
    };
    let eq11 = move |c: &[f64; 15], _c_h2o2: &f64, c_ho2min: &f64, _c_scn: &f64| {
        -k_edta * c_ho2min.powf(2.0) * c[2] * c[11]
    };
    let eq12 = move |c: &[f64; 15], c_h2o2: &f64, c_ho2min: &f64, _c_scn: &f64| {
        k_edta * c_ho2min.powf(2.0) * c[2] * c[11] - k_oh2 * c_h2o2 * c[12]
            + k_edta2 * c_ho2min.powf(2.0) * c[13]
            - k_prim_oxid * c[9] * c[12]
    };
    let eq13 = move |c: &[f64; 15], c_h2o2: &f64, c_ho2min: &f64, _c_scn: &f64| {
        k_oh2 * c_h2o2 * c[12] - k_edta2 * c_ho2min.powf(2.0) * c[13] - k_prim_oxid * c[9] * c[13]
    };
    let eq14 = move |c: &[f64; 15], c_h2o2: &f64, _c_ho2min: &f64, _c_scn: &f64| {
        -k_phen1 * c[5] * c[5] * c[14] - k_phen2 * c_h2o2 * c[14]
    };
    let mut eqs = Vec::new();

    pushinator(|_, _, _, _| 0.0, &mut eqs); // niczemu nie służy, tylko po to,
                                            // żeby liczba elementów była taka sama jak w conc poniżej
    pushinator(eq1, &mut eqs);
    pushinator(eq2, &mut eqs);
    pushinator(eq3, &mut eqs);
    pushinator(eq4, &mut eqs);
    pushinator(eq5, &mut eqs);
    pushinator(eq6, &mut eqs);
    pushinator(eq7, &mut eqs);
    pushinator(eq8, &mut eqs);
    pushinator(eq9, &mut eqs);
    pushinator(eq10, &mut eqs);
    pushinator(eq11, &mut eqs);
    pushinator(eq12, &mut eqs);
    pushinator(eq13, &mut eqs);
    pushinator(eq14, &mut eqs);

    //stezenia
    let mut conc = [1e-8; 15];
    conc[14] = 0.0;
    conc[13] = 0.0;
    conc[12] = 0.0;
    conc[11] = 0.0;
    conc[4] = 8e-5;
    let c_h2o2 = 0.3425;
    let c_ho2min = 0.0375;
    let c_scn = 0.025;

    let mut pot = potencjal_mieszany(c_ho2min, conc[5], conc[2], conc[1]);
    //println!("{t},{},{}", (pot.0), (pot.1));

    //stale do rk4
    let mut k1s = [0.0; 15];
    let mut k2concs = [0.0; 15];
    let mut k2s = [0.0; 15];
    let mut k3concs = [0.0; 15];
    let mut k3s = [0.0; 15];
    let mut k4concs = [0.0; 15];
    let mut k4s = [0.0; 15];

    let f = OpenOptions::new()
        .append(true)
        .open(path1win)
        .expect("Unable to open file");
    let mut f = BufWriter::new(f);
    let stezenia_plik = OpenOptions::new()
        .append(true)
        .open(path2win)
        .expect("Unable to open file");
    let mut stezenia_plik = BufWriter::new(stezenia_plik);
    f.write_all(format!("{},{},{}\n", t / 60.0, pot.0, pot.1).as_bytes())
        .expect("tragedia");

    let zapisy_na_sekunde = 1.0;
    let mut switch = true;
    while t < 10000.0 {
        t += h;
        if t > 2400.0 && switch {
            h = 1e-9;
            conc[14] = 3000e-5;
            switch = false;
        }
        rk4(
            &mut conc,
            &mut h,
            &eqs,
            c_h2o2,
            c_ho2min,
            c_scn,
            &mut k1s,
            &mut k2concs,
            &mut k2s,
            &mut k3concs,
            &mut k3s,
            &mut k4concs,
            &mut k4s,
        );
        pot = potencjal_mieszany(c_ho2min, conc[5], conc[2], conc[1]);

        if (zapisy_na_sekunde * (t + h)).floor() >= (zapisy_na_sekunde * t).ceil() {
            //println!("{t},{},{},{:?}", (pot.0), (pot.1), conc);
            f.write_all(format!("{},{},{}\n", t / 60.0, pot.0, pot.1).as_bytes())
                .expect("tragedia");
            stezenia_plik
                .write_all(
                    format!(
                        "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                        t / 60.0,
                        conc[1],
                        conc[2],
                        conc[3],
                        conc[4],
                        conc[5],
                        conc[6],
                        conc[7],
                        conc[8],
                        conc[9],
                        conc[10],
                        conc[11],
                        conc[12],
                        conc[13],
                        conc[14]
                    )
                    .as_bytes(),
                )
                .expect("tragedia stezenia");
        }
    }
}

fn rk4(
    conc: &mut [f64; 15],
    h: &mut f64,
    eqs: &[BoxedFunction],
    c_h2o2: f64,
    c_ho2min: f64,
    c_scn: f64,
    k1s: &mut [f64; 15],
    k2concs: &mut [f64; 15],
    k2s: &mut [f64; 15],
    k3concs: &mut [f64; 15],
    k3s: &mut [f64; 15],
    k4concs: &mut [f64; 15],
    k4s: &mut [f64; 15],
) {
    let kxconculator = |kxs: &[f64; 15], multiplier, kxconcs: &mut [f64; 15]| {
        (0usize..15)
            .zip(kxs)
            .for_each(|(i, k)| kxconcs[i] = conc[i] + k * multiplier);
    };
    let kxer = |kxconcs: &[f64; 15], kxs: &mut [f64; 15]| {
        (0usize..15).for_each(|i| {
            kxs[i] = *h * (eqs[i].f)(kxconcs, &c_h2o2, &c_ho2min, &c_scn);
        });
    };
    kxer(conc, k1s);
    kxconculator(k1s, 0.5, k2concs);
    kxer(k2concs, k2s);
    kxconculator(k2s, 0.5, k3concs);
    kxer(k3concs, k3s);
    kxconculator(k3s, 1.0, k4concs);
    kxer(k4concs, k4s);
    let rk3_of_c9 = conc[9] + (k1s[9] + 4.0 * k2s[9] + k3s[9]) / 6.0;
    let delta = (conc[9] - rk3_of_c9).abs();
    let mut new_h = *h;
    if delta > 10.0_f64.powf(-7.0) && *h > 10.0_f64.powf(-8.0) {
        new_h = *h / 2.0
    } else if delta < 10.0_f64.powf(-10.0) && *h < 10.0_f64.powf(-3.0) {
        new_h = 2.0 * *h
    }
    if *h != new_h {
        println!("{}", new_h);
    }
    *h = new_h;
    zip!(k1s, k2s, k3s, k4s)
        .enumerate()
        .for_each(|(id, (k1, (k2, (k3, k4))))| {
            conc[id] += (k1 + 2.0 * k2 + 2.0 * k3 + *k4) / 6.0;
        });
}

fn potencjal_mieszany(c_ho2min: f64, c_ho2rod: f64, c_cuoh3: f64, c_cuoh2: f64) -> (f64, f64) {
    const F: f64 = 96485.3321;
    const R: f64 = 8.314462;
    const T: f64 = 283.15;

    let e1 = -0.18 + R * T / F * (c_cuoh3 / c_cuoh2).ln();
    let e2 = 0.22 + R * T / F * (c_ho2rod / c_ho2min).ln();
    let i01 = F * c_cuoh3.powf(0.5) * c_cuoh2.powf(0.5);
    let i02pt = 10.0_f64.powf(-5.0) * F * c_ho2min.powf(0.5) * c_ho2rod.powf(0.5);
    let i02au = 10.0_f64.powf(-8.0) * F * c_ho2min.powf(0.5) * c_ho2rod.powf(0.5);

    let e_pt = R * T / F
        * ((i01 * (F * e1 / 2.0 / R / T).exp() + i02pt * (F * e2 / 2.0 / R / T).exp())
            / (i01 * (-F * e1 / 2.0 / R / T).exp() + i02pt * (-F * e2 / 2.0 / R / T).exp()))
        .ln();
    let e_au = R * T / F
        * ((i01 * (F * e1 / 2.0 / R / T).exp() + i02au * (F * e2 / 2.0 / R / T).exp())
            / (i01 * (-F * e1 / 2.0 / R / T).exp() + i02au * (-F * e2 / 2.0 / R / T).exp()))
        .ln();
    (e_au, e_pt)
}
