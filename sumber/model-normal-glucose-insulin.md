# Model Glukosa-Insulin: Mode Normal (Pankreas Sehat)

Dokumentasi ini merangkum model matematis, sumber, parameter, modifikasi, dan hasil
verifikasi untuk **Mode Normal** pada simulator sistem glukosa-insulin. Dokumen ini
menjadi fondasi sebelum pengembangan Mode Diabetes (tanpa kontrol dan dengan
Artificial Pancreas/PID) serta visualisasi interaktif.

> **Batas penggunaan:** simulator ini dibuat hanya untuk pendidikan dan
> visualisasi konsep dinamika glukosa-insulin. Seluruh kurva, dosis, dan skenario
> merupakan keluaran model matematis populasi, bukan prediksi klinis individual,
> alat diagnosis, rekomendasi terapi, atau dasar penentuan dosis insulin bagi
> pasien nyata.

---

## 1. Ringkasan Model

Mode Normal merepresentasikan kondisi tubuh sehat, di mana pankreas secara alami
mensekresi insulin sebagai respons terhadap kenaikan glukosa darah. Sistem ini
merupakan **closed-loop** (loop tertutup) — glukosa memengaruhi sekresi insulin, dan
insulin balik memengaruhi glukosa.

### Pemetaan ke Block Diagram Kontrol

| Blok | Isi di Mode Normal |
|---|---|
| Input / Setpoint | `p5` (ambang glukosa) |
| Controller | Persamaan `dI/dt` (fungsi sekresi pankreas alami) |
| Process (Plant) | Persamaan `dG/dt` dan `dX/dt` (respons tubuh terhadap insulin) |
| Feedback path | `G(t)` dibaca balik ke suku `[G-p5]⁺` |
| Disturbance | `D(t)` (input makan) |

---

## 2. Persamaan Matematis

Sistem 3 persamaan diferensial nonlinear, closed-loop:

```
dG(t)/dt = -p1(G(t) - Gb) - X(t)·G(t) + D(t)

dX(t)/dt = -p2·X(t) + p3(I(t) - Ib)

dI(t)/dt = p6·[G(t) - p5]⁺ · τ(t) - n(I(t) - Ib)
```

**Keterangan variabel:**
- `G(t)` — konsentrasi glukosa plasma (mg/dL)
- `X(t)` — efek insulin aktif di jaringan remote (1/menit)
- `I(t)` — konsentrasi insulin plasma (μU/mL)
- `D(t)` — input gangguan makan (rate of glucose appearance)
- `τ(t)` — waktu lokal sejak G(t) terakhir melewati ambang p5 (lihat Bagian 4)
- `[x]⁺` — operator `max(x, 0)`

**Kondisi awal (baseline/istirahat):** `G(0)=Gb`, `X(0)=0`, `I(0)=Ib`

---

## 3. Sumber (Rantai Sitasi)

| Bagian Model | Sumber |
|---|---|
| Persamaan glukosa (`dG/dt`, `dX/dt`) — "Model VI" | Bergman, R.N., Ider, Y.Z., Bowden, C.R., & Cobelli, C. (1979). Quantitative estimation of insulin sensitivity. *American Journal of Physiology*, 236(6), E667–E677. |
| Persamaan insulin (`dI/dt`) + nilai parameter contoh subjek normal | Pacini, G., & Bergman, R.N. (1986). MINMOD: a computer program to calculate insulin sensitivity and pancreatic responsivity from the frequently sampled intravenous glucose tolerance test. *Computer Methods and Programs in Biomedicine*, 23(2), 113–122. |
| Penggabungan (*coupling*) kedua model jadi sistem closed-loop 3-persamaan; identifikasi bug equilibrium | Friis-Jensen, E. (2007). *Modeling and Simulation of Glucose-Insulin Metabolism*. B.Sc. Thesis, Technical University of Denmark (DTU). |
| Fungsi input makan `D(t) = A·t·e^(-t/τ)/τ²` (absorpsi 2-kompartemen berantai, konstanta transfer identik) | Hovorka, R., Canonico, V., Chassin, L.J., Haueter, U., Massi-Benedetti, M., Orsini Federici, M., Pieber, T.R., Schaller, H.C., Schaupp, L., Vering, T., & Wilinska, M.E. (2004). Nonlinear model predictive control of glucose concentration in subjects with type 1 diabetes. *Physiological Measurement*, 25(4), 905–920 — Persamaan (4). Nilai τ = tmax,G = 40 menit (Tabel 1 paper yang sama). |

**Catatan penting:** Persamaan glukosa dan insulin *aslinya* dirancang **open-loop**
untuk analisis data IVGTT (Intravenous Glucose Tolerance Test) — insulin/glukosa
terukur dipakai sebagai *input*, bukan hasil simulasi timbal-balik. Coupling menjadi
sistem closed-loop adalah langkah tambahan (Friis-Jensen, 2007), bukan bagian dari
paper asli Bergman (1979) maupun Pacini & Bergman (1986).

---

## 4. Parameter

### Parameter Model (dari Pacini & Bergman, 1986 — Fig. 3A/3B, dataset "NORMAL.DAT")

| Simbol | Nilai | Satuan | Arti Fisiologis |
|---|---|---|---|
| p1 | 0.03082 | 1/min | Glucose effectiveness (SG) |
| p2 | 0.02093 | 1/min | Peluruhan efek insulin aktif |
| p3 | 1.062×10⁻⁵ | L/(min²·mU) | Insulin sensitivity |
| n | 0.3 | 1/min | Clearance insulin plasma |
| p6 (γ) | 0.003349 | — | Laju respons pankreas |
| Gb | 92 | mg/dL | Glukosa basal |
| Ib | 7.3 | μU/mL | Insulin basal |

### Parameter yang DIMODIFIKASI dari nilai asli paper

| Simbol | Nilai Asli Paper | Nilai Dipakai Simulator | Alasan |
|---|---|---|---|
| p5 | 89.5 | **94** | p5 asli < Gb menyebabkan sistem tidak bisa mencapai equilibrium (lihat Bagian 5) |
| Suku waktu di dI/dt | `t` (absolut, sejak simulasi/eksperimen mulai) | **τ (lokal, reset tiap episode)** | `t` absolut menyebabkan sekresi insulin tumbuh tak terbatas / tidak konsisten pada simulasi jangka panjang dengan banyak episode (lihat Bagian 5) |

---

## 5. Modifikasi & Justifikasi

### 5.1 Modifikasi p5: 89.5 → 94

**Masalah:** Pembuktian analitik menunjukkan bahwa titik equilibrium (G=Gb, X=0, I=Ib)
hanya valid jika `p5 ≥ Gb`:

```
Di titik (G=Gb, I=Ib):  dI/dt = p6·[Gb - p5]⁺·τ

Jika p5 ≥ Gb  →  [Gb-p5]⁺ = 0  →  dI/dt = 0   ✓ equilibrium tercapai
Jika p5 < Gb  →  [Gb-p5]⁺ > 0  →  dI/dt ≠ 0   ✗ equilibrium TIDAK tercapai
```

Nilai asli paper (p5=89.5) lebih kecil dari Gb=92, sehingga sistem tidak akan pernah
diam meski tanpa gangguan apapun — insulin akan terus "bocor" naik. Ini tidak masalah
di paper asli (yang hanya fit data 182 menit sekali jalan), tapi jadi masalah untuk
simulasi closed-loop kontinu.

**Solusi:** p5 dinaikkan menjadi 94 (> Gb=92).

### 5.2 Modifikasi suku waktu: t absolut → τ lokal

**Masalah:** Suku `p6·[G-p5]⁺·t` merepresentasikan sekresi insulin fase-kedua
(second-phase), yang seharusnya bergantung pada **berapa lama episode glukosa
tinggi SAAT INI berlangsung** — bukan "berapa lama sejak simulasi dimulai". Jika
memakai `t` absolut, makan yang identik namun terjadi di menit ke-10 vs menit
ke-800 dalam simulasi menghasilkan respons insulin yang jauh berbeda (dibuktikan:
29,5 vs 102,3 μU/mL — beda 3,5×), padahal gangguannya sama persis.

**Solusi:** Ganti `t` dengan `τ` — waktu lokal yang **direset ke 0** setiap kali
`G(t)` turun ≤ p5, dan bertambah selama `G(t)` > p5.

---

## 6. Verifikasi

### 6.1 Verifikasi Analitik (Syarat Equilibrium)
Dibuktikan aljabar bahwa `p5 ≥ Gb` adalah syarat perlu dan cukup (bersama `dG/dt`
dan `dX/dt` yang otomatis nol di titik basal) agar (Gb, 0, Ib) menjadi titik
equilibrium yang valid.

### 6.2 Verifikasi Numerik — Uji Diam (Lapis 2)
Simulasi 1000 menit dari kondisi awal G(0)=Gb, X(0)=0, I(0)=Ib, **tanpa gangguan**.

| Variabel | Deviasi Maksimum dari Baseline |
|---|---|
| G(t) | 0.000000 mg/dL |
| X(t) | 0.00000000 |
| I(t) | 0.000000 μU/mL |

✅ **Lolos** — sistem benar-benar diam, tidak drift.

### 6.3 Verifikasi Numerik — Uji Gangguan lalu Kembali (Lapis 3)
Simulasi 1000 menit dengan 1× trigger makan di menit ke-60.

| Metrik | Hasil |
|---|---|
| Glukosa puncak | 122.55 mg/dL @ t=120.5 menit |
| Insulin puncak | 27.13 μU/mL @ t=142.0 menit (delay 21.5 menit dari puncak glukosa) |
| Glukosa akhir (t=1000 menit) | 92.0000 mg/dL (selisih dari Gb: 0.000002 mg/dL) |

✅ **Lolos** — kembali stabil sempurna ke baseline, termasuk pola *undershoot*
realistis sebelum stabil (konsisten dengan fenomena yang dilaporkan di Bergman
et al., 1979).

### 6.4 Verifikasi Konsistensi τ Lokal (Time-Shift Invariance)
Makan identik ditrigger di 6 titik waktu berbeda (t=10, 100, 300, 500, 800, 1500
menit).

| Titik Waktu Makan | G Puncak | I Puncak |
|---|---|---|
| Semua titik (10–1500) | 122.55 mg/dL | 27.13 μU/mL |

✅ **Lolos** — hasil identik di semua titik waktu, membuktikan τ lokal
menghilangkan ketergantungan artifisial pada "jam berapa sekarang".

### 6.5 Validasi Implementasi terhadap Data Paper
Model glukosa (`dG/dt`, `dX/dt`) diuji **open-loop** (insulin = data pengukuran asli
dari paper, bukan hasil simulasi kita) dan dibandingkan dengan hasil fit resmi paper
(Fig. 3A, kolom COMP.GLU).

**RMSE = 1.36 mg/dL** — sangat kecil, mengonfirmasi implementasi model glukosa
kita benar dan sesuai paper.

Saat diuji **closed-loop penuh** (insulin dari simulasi kita sendiri via τ lokal)
dibandingkan data mentah (dengan alignment waktu +4 menit dan eksklusi 4 titik
awal yang di-zero-weight oleh paper karena fase mixing ekstraseluler):

**RMSE = 10.02 mg/dL** — selisih ini wajar dan diharapkan, karena fungsi sekresi
pankreas adalah model/pendekatan matematis, bukan replikasi data insulin individu
pasien tertentu.

### Ringkasan Status Verifikasi

| # | Uji | Hasil |
|---|---|---|
| 1 | Analitik (syarat equilibrium) | ✅ Terbukti aljabar |
| 2 | Numerik — uji diam | ✅ Deviasi 0.000000 |
| 3 | Numerik — uji kembali stabil | ✅ Kembali ke baseline (selisih 0.000002) |
| 4 | Konsistensi τ lokal | ✅ Identik di semua titik waktu |
| 5 | Validasi vs data paper (model glukosa) | ✅ RMSE 1.36 mg/dL |
| 6 | Validasi vs data paper (closed-loop penuh) | ⚠️ RMSE 10.02 mg/dL (wajar) |

---

## 7. Perilaku Fisiologis yang Berhasil Direproduksi

- Glukosa naik ~30 mg/dL setelah makan, puncak ~60 menit pasca-makan
- Insulin menyusul naik dengan delay ~21.5 menit dari puncak glukosa (konsisten
  dengan pola first/second-phase insulin release)
- Sedikit *undershoot* sebelum stabil kembali ke baseline
- Kembali stabil sempurna tanpa osilasi liar, baik setelah gangguan tunggal maupun
  ketika diuji pada berbagai titik waktu simulasi

---

## 8. Verifikasi Sumber Fungsi Input Makan D(t)

Fungsi `D(t)` diverifikasi langsung ke PDF paper asli. Persamaan (4) di
Hovorka dkk. (2004):

```
UG(t) = DG · AG · t · e^(-t/tmax,G) / t²max,G
```

Perbandingan suku-per-suku dengan kode simulator:

| Hovorka (2004) | Simulator | Status |
|---|---|---|
| `DG·AG` (karbohidrat × bioavailabilitas) | `A` (amplitudo porsi makan) | Setara secara bentuk, digabung menjadi satu konstanta |
| `t` | `dtSince` (waktu sejak makan) | Identik |
| `e^(-t/tmax,G)` | `e^(-dtSince/τ_meal)` | Identik |
| `1/t²max,G` (normalisasi) | `1/τ_meal²` | Identik |
| `tmax,G = 40 min` | `mealTau = 40` | Identik |

Bentuk tersebut adalah solusi analitik dua kompartemen berantai dengan laju
transfer identik:

```
dD1/dt = cho - D1/tmax
dD2/dt = D1/tmax - D2/tmax
UG     = D2/tmax
```

**Catatan koreksi sumber:**

- Fisher (1991) tidak memakai bentuk ini; paper tersebut memakai
  `P(t) = B·exp(-kt)`, yaitu eksponensial sederhana.
- Urbina dkk. (2020) memakai `a0·exp(-b0·t)`.
- Lehmann & Deutsch (1992) memakai fungsi pengosongan lambung
  trapezoidal/segitiga.
- Hovorka dkk. (2004) adalah sumber yang cocok secara literal dengan bentuk
  fungsi yang digunakan simulator.

**Batas verifikasi:** bentuk fungsi dan `mealTau = 40` menit telah didukung
oleh Hovorka, tetapi amplitudo `A = 80/150/220` masih berupa skala kasar agar
respons terlihat wajar pada gabungan parameter Bergman–Pacini. Nilai tersebut
bukan konversi gram karbohidrat. Hovorka menggunakan `DG·AG` dalam satuan
model yang berbeda.

**Exercise tidak menjadi bagian model.** Perubahan sementara `p3` yang pernah
direncanakan belum memiliki rujukan dan validasi yang memadai, sehingga tidak
boleh disebut sebagai representasi fisiologis olahraga. Exercise hanya dapat
ditambahkan kembali setelah memiliki model, parameter, dan pengujian tersendiri.

---

## 9. Parameter Final untuk Implementasi Visualisasi

```python
p1 = 0.03082
p2 = 0.02093
p3 = 1.062e-5
n  = 0.3
p6 = 0.003349
Gb = 92.0
Ib = 7.3
p5 = 94.0        # dimodifikasi dari 89.5
G0 = 92.0        # mulai dari baseline (bukan G0=287 ala IVGTT)
X0 = 0.0
I0 = 7.3         # = Ib, mulai dari baseline

# Suku waktu pakai tau lokal (reset saat G <= p5), BUKAN t absolut
```

---

## 10. Rancangan Mode Diabetes Tipe 1 Tanpa PID

> **Status:** spesifikasi T1D v2 telah melewati verifikasi matematis dan numerik
> terpisah, tetapi belum diimplementasikan pada simulator UI. Mode ini sengaja
> dibatasi pada terapi insulin manual (*open-loop*), bukan Artificial Pancreas
> atau PID.

### 10.1 Tujuan dan ruang lingkup

Mode ini merepresentasikan penderita diabetes tipe 1 (T1D) yang tidak memiliki
sekresi insulin endogen yang memadai. Insulin diberikan dari luar tubuh sebagai:

- **basal**, yaitu infus konstan untuk kebutuhan insulin di antara waktu makan;
- **bolus makan**, yaitu dosis manual pada waktu yang dipilih pengguna; dan
- **bolus koreksi**, yaitu dosis manual untuk mengoreksi hiperglikemia.

Tidak ada hubungan otomatis dari pengukuran glukosa ke dosis insulin. Karena itu,
mode ini disebut *open-loop*. Pengguna dapat memilih makan tanpa bolus, memberi
bolus tepat waktu, memberi bolus terlambat, atau memberi dosis berlebih.

Mode ini bukan model ketoasidosis diabetik dan tidak ditujukan untuk menghitung
dosis insulin klinis bagi pasien nyata.

### 10.2 Dasar ilmiah terapi basal-bolus

Terapi penggantian insulin pada T1D secara fisiologis terdiri dari insulin basal,
insulin waktu makan, dan insulin koreksi. Terapi ini diberikan melalui beberapa
suntikan harian atau *continuous subcutaneous insulin infusion* (CSII), yaitu
pompa yang memberikan basal secara kontinu dan bolus makan secara manual.

Dasar peer-reviewed:

1. Holt, R.I.G., DeVries, J.H., Hess-Fischl, A., et al. (2021). *The management
   of type 1 diabetes in adults: A consensus report by the American Diabetes
   Association (ADA) and the European Association for the Study of Diabetes
   (EASD).* Diabetologia, 64, 2609-2652.
   DOI: https://doi.org/10.1007/s00125-021-05568-3
2. The Diabetes Control and Complications Trial Research Group (1993). *The
   effect of intensive treatment of diabetes on the development and progression
   of long-term complications in insulin-dependent diabetes mellitus.* New
   England Journal of Medicine, 329, 977-986.
   DOI: https://doi.org/10.1056/NEJM199309303291401

Kedua sumber tersebut mendukung keberadaan input insulin manual pada simulator.
Nilai dosis pada simulator tetap harus dikalibrasi terhadap model dan tidak boleh
ditafsirkan sebagai rekomendasi medis.

### 10.3 Keputusan struktur model T1D v2

Percobaan awal T1D v1 memakai model glukosa `G-X` Mode Normal dan hanya mengganti
sekresi pankreas dengan absorpsi insulin subkutan. Model lima-state tersebut
stabil secara matematis, tetapi gagal membedakan keadaan T1D dari Mode Normal:

- penghentian insulin hanya menghasilkan equilibrium sekitar `104.57 mg/dL`; dan
- makan `A = 150` tanpa bolus hanya mencapai sekitar `127.14 mg/dL`, sedangkan
  Mode Normal mencapai sekitar `122.55 mg/dL`.

Penyebabnya adalah suku `-p1(G-Gb)`, yang tetap menarik glukosa kembali ke
baseline walaupun insulin hampir tidak tersedia. Karena itu, T1D v2 tidak lagi
memakai persamaan glukosa Bergman-Pacini Mode Normal. T1D v2 memakai seluruh
subsistem glukosa, insulin, dan aksi insulin dari Hovorka dkk. (2004).

Mode Normal tetap memakai persamaannya sendiri dan tidak diubah oleh keputusan
ini.

### 10.4 Konvensi state dan satuan

T1D v2 memiliki delapan state:

```text
Q1, Q2, S1, S2, I, x1, x2, x3
```

| State | Arti |
|---|---|
| `Q1` | Massa glukosa pada kompartemen plasma/accessible |
| `Q2` | Massa glukosa pada kompartemen jaringan/non-accessible |
| `S1`, `S2` | Depot absorpsi insulin subkutan |
| `I` | Konsentrasi insulin plasma |
| `x1` | Efek insulin pada distribusi/transport glukosa |
| `x2` | Efek insulin pada disposal glukosa |
| `x3` | Efek insulin pada produksi glukosa endogen |

Implementasi awal sebaiknya memakai bentuk ternormalisasi per kilogram seperti
paper: `Q1` dan `Q2` dalam mmol/kg, `VG` dan `VI` dalam L/kg, laju glukosa dalam
mmol/kg/menit, insulin dalam mU/L, serta laju insulin dalam mU/kg/menit.
`1 mU/L` bernilai numerik sama dengan `1 microU/mL`.

Konsentrasi glukosa internal dihitung dalam mmol/L:

```text
G(t) = Q1(t) / VG
```

Untuk tampilan UI:

```text
G_mg_dL(t) = 18.0182 * G_mmol_L(t)
```

### 10.5 Subsistem glukosa Hovorka

```text
dQ1/dt = -F01c - x1*Q1
          + k12*Q2 - FR + UG + EGP0*max(0, 1 - x3)

dQ2/dt = x1*Q1 - (k12 + x2)Q2

G = Q1/VG
```

Paper Hovorka menulis suku pertama dalam bentuk literal
`-[F01c/(VG*G) + x1]Q1`. Karena `Q1 = VG*G`, bagian
`-[F01c/(VG*G)]Q1` identik dengan `-F01c` ketika `G > 0`. Implementasi memakai
bentuk langsung `-F01c` agar tidak menghasilkan `0/0` saat `G = 0`.

`F01c` adalah penggunaan glukosa yang tidak bergantung pada insulin dan dikoreksi
pada konsentrasi glukosa rendah:

```text
F01c = F01                  jika G >= 4.5 mmol/L
F01c = F01*(G/4.5)          jika G <  4.5 mmol/L
```

`FR` adalah ekskresi glukosa melalui ginjal, yang aktif di atas ambang
`9 mmol/L` atau sekitar `162 mg/dL`:

```text
FR = 0.003*(G - 9)*VG       jika G >= 9 mmol/L
FR = 0                      jika G <  9 mmol/L
```

Suku produksi glukosa endogen pada paper adalah `EGP0*(1-x3)`. Bentuk literal
tersebut dapat menjadi negatif ketika `x3 > 1`, sehingga pada bolus ekstrem hati
secara numerik berubah menjadi penyerap glukosa dan massa `Q1` dapat melewati
nol. Proyek ini menerapkan batas:

```text
EGP = EGP0*max(0, 1 - x3)
```

Clipping tersebut adalah keputusan implementasi proyek untuk menjaga produksi
glukosa endogen tetap non-negatif; bukan bagian yang tertulis pada persamaan
Hovorka (2004). Pengaman ini tidak menambahkan hormon kontraregulasi dan tidak
mencegah hipoglikemia berat akibat overdosis insulin.

### 10.6 Input makanan: bentuk sama, skala dan satuan berbeda

Input makanan T1D menggunakan persamaan (4) Hovorka:

```text
tauMeal = t - tMeal

UG(t) = DG*AG*tauMeal*exp(-tauMeal/tmaxG)/(tmaxG^2), jika tauMeal >= 0
UG(t) = 0,                                             jika tauMeal < 0
```

Persamaan tersebut mempunyai bentuk kernel yang sama dengan fungsi Mode Normal:

```text
D(t) = A*tauMeal*exp(-tauMeal/tau)/(tau^2)
```

Pemetaan bentuknya adalah:

```text
A   <-> DG*AG
tau <-> tmaxG
```

Jadi, bentuk kurva makanan lama memang berasal dari Hovorka. Perbedaannya bukan
pada bentuk persamaan, tetapi pada definisi amplitudo, satuan, dan state tujuan:

- Mode Normal memakai `A = 80/150/220` sebagai amplitudo internal dan `D(t)`
  langsung masuk ke `dG/dt`;
- Mode T1D memakai `DG` sebagai dosis karbohidrat dalam mmol/kg, mengalikan
  bioavailabilitas `AG`, dan memasukkan `UG` ke persamaan massa `dQ1/dt`.

Konversi gram karbohidrat menjadi `DG`, dengan pendekatan karbohidrat sebagai
ekuivalen glukosa, adalah:

```text
DG [mmol/kg] = carbsGram*1000/(180.16*weightKg)
```

Parameter Hovorka:

```text
AG    = 0.8
tmaxG = 40 min
```

Untuk beberapa kejadian makan, `UG_total(t)` adalah jumlah kontribusi semua
makanan yang sudah dimulai. Input `80/150/220` milik Mode Normal tidak boleh
dianggap sebagai gram karbohidrat dan tidak boleh dipakai langsung sebagai `DG`.

### 10.7 Subsistem absorpsi dan eliminasi insulin

```text
dS1/dt = u(t) - S1/tmaxI

dS2/dt = S1/tmaxI - S2/tmaxI

UI(t)  = S2/tmaxI

dI/dt  = UI/VI - ke*I
```

Input insulin manual:

```text
u(t) = uBasal(t) + uBolus(t)
```

`uBasal` bersifat kontinu. `uBolus` adalah dosis selama interval pendek pada
waktu yang dipilih pengguna. Makan tidak otomatis memberikan bolus. Tidak ada
sekresi insulin endogen dan tidak ada jalur `G -> pankreas -> I`.

Jika basal konstan, equilibrium subsistem insulin memenuhi:

```text
S1b = uBasal*tmaxI
S2b = uBasal*tmaxI
Ib  = uBasal/(VI*ke)
```

Untuk bolus `B` unit pada subjek dengan berat `W` kg, jumlah insulin yang masuk
ke depot adalah:

```text
bolus_mU_per_kg = B*1000/W
```

Jika bolus diterapkan selama satu langkah waktu `dt`, laju inputnya harus dibagi
dengan `dt` agar integral dosis tetap sama saat `dt` berubah.

### 10.8 Subsistem aksi insulin

```text
dx1/dt = -ka1*x1 + kb1*I
dx2/dt = -ka2*x2 + kb2*I
dx3/dt = -ka3*x3 + kb3*I
```

Paper memberikan sensitivitas alternatif `SIT`, `SID`, dan `SIE`, dengan:

```text
kb1 = ka1*SIT
kb2 = ka2*SID
kb3 = ka3*SIE
```

`x1`, `x2`, dan `x3` harus bereaksi terhadap insulin plasma yang telah melalui
penundaan absorpsi `S1-S2`; bolus tidak boleh memengaruhi glukosa secara instan.

### 10.9 Parameter populasi awal

| Parameter | Nilai awal | Satuan |
|---|---:|---|
| `k12` | 0.066 | /menit |
| `ka1` | 0.006 | /menit |
| `ka2` | 0.060 | /menit |
| `ka3` | 0.030 | /menit |
| `ke` | 0.138 | /menit |
| `VI` | 0.12 | L/kg |
| `VG` | 0.16 | L/kg |
| `AG` | 0.8 | tanpa satuan |
| `tmaxG` | 40 | menit |
| `tmaxI` | 55 | menit |
| `SIT` | `51.2e-4` | /menit per mU/L |
| `SID` | `8.2e-4` | /menit per mU/L |
| `SIE` | `520e-4` | per mU/L |
| `EGP0` | 0.0161 | mmol/kg/menit |
| `F01` | 0.0097 | mmol/kg/menit |

Nilai tersebut adalah nilai populasi awal dari Hovorka, bukan dosis personal
untuk pengguna. Berat badan diperlukan untuk mengonversi input gram dan unit
menjadi besaran per kilogram.

Sumber model dan parameter:

- Hovorka, R., Canonico, V., Chassin, L.J., et al. (2004). *Nonlinear model
  predictive control of glucose concentration in subjects with type 1 diabetes.*
  Physiological Measurement, 25(4), 905-920.
  DOI: https://doi.org/10.1088/0967-3334/25/4/010

### 10.10 Kondisi awal dan equilibrium

Simulasi tidak boleh dimulai dengan menebak `Q1`, `Q2`, atau state aksi insulin
secara terpisah. Untuk target glukosa puasa dan basal yang dipilih:

1. hitung equilibrium `S1`, `S2`, dan `I` dari basal;
2. hitung equilibrium `x1`, `x2`, dan `x3` dari `dxi/dt = 0`;
3. selesaikan `dQ1/dt = 0` dan `dQ2/dt = 0` tanpa makanan untuk `Q1` dan `Q2`;
4. pastikan `Q1/VG` sama dengan target glukosa yang konsisten dengan basal; dan
5. gunakan hasil equilibrium sebagai seluruh state awal.

Jika target glukosa dan basal tidak menghasilkan equilibrium yang sama, salah
satunya harus dicari secara numerik. Memaksa keduanya secara terpisah akan
menciptakan transien palsu pada awal simulasi.

### 10.11 Skenario demonstrasi yang dituju

| Skenario | Input pengguna | Perilaku yang diharapkan |
|---|---|---|
| Basal stabil | Basal aktif, tanpa makan | Sistem tetap di sekitar equilibrium |
| Makan tanpa bolus | Karbohidrat, bolus 0 | Glukosa meningkat jelas dan bertahan lebih lama |
| Bolus tepat waktu | 50 g karbohidrat dan bolus awal 2 U | Lonjakan glukosa lebih terkendali |
| Bolus terlambat | Makan lebih dahulu | Puncak dan paparan hiperglikemia lebih besar |
| Bolus berlebih | Makanan kecil, bolus besar | Risiko hipoglikemia akibat insulin aktif |
| Basal dihentikan | `uBasal = 0` | Insulin dan aksi insulin turun; produksi glukosa hati meningkat |

Mode ini tidak memodelkan keton, ketoasidosis, dehidrasi, perubahan hormon
kontraregulasi, atau perjalanan penyakit tanpa terapi selama beberapa hari.

Untuk preset demonstrasi awal pada asumsi berat `70 kg`, gunakan bolus `2 U`
untuk makanan `50 g`. Nilai ini dipilih dari perilaku model populasi yang telah
diuji, bukan rekomendasi dosis klinis. Bolus `5 U` tidak boleh menjadi default
karena pada konfigurasi ini menghasilkan hipoglikemia berat.

### 10.12 Verifikasi wajib sebelum implementasi UI

T1D v2 belum dianggap siap untuk UI sebelum melewati seluruh uji berikut:

1. **Equilibrium basal:** seluruh delapan state tidak drift tanpa makanan.
2. **Neraca makanan:** integral `UG` harus sama dengan `DG*AG` dalam toleransi
   numerik.
3. **Neraca bolus:** integral `uBolus` harus sama dengan dosis yang diminta dan
   tidak bergantung pada `dt`.
4. **Penghentian basal:** `I`, `x1`, `x2`, dan `x3` menurun, sedangkan glukosa
   meningkat karena penekanan `EGP` berkurang.
5. **Makan tanpa bolus:** menghasilkan hiperglikemia yang jelas, bukan respons
   yang hampir sama dengan Mode Normal.
6. **Dosis-respons:** bolus yang meningkat menurunkan puncak atau AUC glukosa
   secara konsisten dalam rentang uji.
7. **Waktu bolus:** bolus terlambat menghasilkan paparan hiperglikemia lebih
   besar daripada bolus tepat waktu.
8. **Overdosis:** dosis berlebih dapat menimbulkan hipoglikemia tanpa membuat
   massa glukosa atau insulin negatif.
9. **Fungsi potongan:** uji khusus di sekitar `G = 4.5` dan `G = 9 mmol/L` untuk
   memastikan `F01c` dan `FR` tidak salah cabang.
10. **Konvergensi numerik:** hasil stabil saat `dt` diperkecil.
11. **Uji satuan:** gram, mmol, kg, U, mU, liter, menit, mmol/L, dan mg/dL
    diperiksa secara eksplisit.

### 10.13 Hasil verifikasi numerik T1D v2

Verifikasi dijalankan dengan asumsi:

```text
Berat badan  = 70 kg
Glukosa awal = 5.5 mmol/L
Makanan      = 50 g karbohidrat
```

Dari 13 pemeriksaan yang berasal dari 11 kelompok uji pada Bagian 10.12,
12 pemeriksaan lulus. Satu kegagalan pada persamaan literal menghasilkan
perbaikan implementasi `max(0, 1-x3)` yang didokumentasikan pada Bagian 10.5.

| Pemeriksaan | Hasil | Status |
|---|---|---|
| Equilibrium basal | Delapan state tidak drift; seluruh eigenvalue negatif; konstanta waktu paling lambat sekitar 733 menit | Lolos |
| Neraca makanan | Integral `UG = DG*AG = 3.1718 mmol/kg`; bentuk tertutup dan rantai dua kompartemen cocok | Lolos |
| Neraca bolus | `5 U = 71.43 mU/kg`; integral `UI` sama dan tidak bergantung pada durasi pulsa | Lolos |
| Basal dihentikan | `I` dan `x1-x3` turun; `G` naik monoton menuju sekitar `402 mg/dL` | Lolos |
| Makan tanpa bolus | Puncak sekitar `283 mg/dL`; berada di atas `180 mg/dL` selama sekitar 500 menit | Lolos |
| Dosis-respons 0-4 U | Puncak dan AUC glukosa turun secara konsisten | Lolos |
| Waktu bolus 2 U | Bolus -15/0/+30/+60 menit menghasilkan puncak sekitar 166/183/220/251 mg/dL | Lolos |
| Overdosis, persamaan literal | Bolus 10 U dengan makan 10 g membuat `x3 = 2.82`, `Q1 = -0.30 mmol/kg`, dan `G = -34 mg/dL` | **Gagal** |
| Overdosis dengan clipping EGP | Massa tidak negatif; nadir `1.43 mg/dL` karena kontraregulasi tidak dimodelkan | Lolos untuk invariansi massa |
| Fungsi potongan | `F01c` dan `FR` kontinu pada 4.5 dan 9 mmol/L; cabang benar | Lolos |
| Konvergensi RK4 | Pada `dt = 1`, galat puncak `0.00137 mg/dL` dan galat nadir `2.24e-5 mg/dL` terhadap RK4 `dt = 0.1` | Lolos |
| Pembeda integrator | Euler `dt = 1` memiliki galat puncak sekitar `0.359 mg/dL` | Lolos |
| Satuan | Seluruh konversi yang diuji konsisten | Lolos |

Hasil tersebut menunjukkan bahwa struktur T1D v2 sudah jauh lebih sesuai untuk
simulator edukasi daripada T1D v1. Penghentian insulin dan makan tanpa bolus kini
memberikan perilaku yang jelas berbeda dari Mode Normal.

#### Catatan kalibrasi preset

Dengan parameter populasi Hovorka dan berat `70 kg`, respons model untuk makanan
`50 g` menjadi terlalu agresif pada bolus yang lebih besar:

- rentang yang diuji masih aman sampai sekitar `2.5 U`;
- `3 U` menghasilkan nadir sekitar `57 mg/dL`; dan
- `5 U` menghasilkan nadir `18.72 mg/dL` dengan persamaan EGP literal atau
  `27.34 mg/dL` setelah clipping EGP.

Karena itu, preset UI direncanakan memakai `2 U`, bukan rasio dosis yang diambil
dari praktik klinis umum. Nilai simulator tidak boleh digunakan untuk menentukan
terapi nyata.

#### Reproduksibilitas kondisi awal

Seluruh state awal wajib berasal dari penyelesaian equilibrium Bagian 10.10.
Fungsi implementasi seperti `equilibrium()` harus menghitung state tersebut dan
tidak menggantinya dengan nilai tebakan.

### 10.14 Status implementasi

- Model hibrida T1D v1 lima-state dinyatakan **tidak digunakan**.
- Spesifikasi T1D v2 Hovorka delapan-state sudah ditetapkan dalam dokumen ini.
- T1D v2 belum diimplementasikan pada `script.js`.
- Hasil verifikasi numerik T1D v2 sudah didokumentasikan pada Bagian 10.13.
- Verifikasi dapat dijalankan ulang dengan `node verification/t1d-v2.js`.
- Skrip `verification/t1d-v1.js` hanya merekam percobaan lama dan bukan
  verifikasi T1D v2.
- Validasi fisiologis terhadap data pasien dapat ditunda untuk tahap berikutnya,
  tetapi verifikasi matematis, satuan, dan perilaku dasar di atas tetap wajib.
- Plant T1D v2 menjadi dasar untuk rancangan Artificial Pancreas pada Bagian 11.

## 11. Rancangan Artificial Pancreas PID-IFB untuk T1D

> **Status:** model PID-IFB telah dituning untuk plant proyek dan melewati
> verifikasi closed-loop numerik pada skenario Bagian 11.13, tetapi belum
> diimplementasikan pada UI. Artificial Pancreas (AP) di bagian ini hanya untuk
> demonstrasi edukasi sistem kontrol, bukan algoritma terapi atau perangkat medis.

### 11.1 Definisi sistem

Sistem disebut Artificial Pancreas karena membentuk loop otomatis:

```text
T1D plant -> glukosa -> sensor CGM -> PID-IFB -> pompa -> insulin -> T1D plant
```

Komponennya adalah:

1. plant T1D v2 Hovorka delapan-state dari Bagian 10;
2. model pembacaan sensor glukosa `Gs`;
3. controller PID diskrit;
4. estimator insulin untuk insulin feedback (IFB); dan
5. batas keluaran pompa dan pengaman integral.

Jika controller tidak menerima jumlah karbohidrat atau bolus makan manual,
skenarionya adalah *fully closed-loop* dalam simulasi. Jika jumlah makan atau
bolus diberikan kepada controller, sistem menjadi *hybrid closed-loop*.

Versi pertama dirancang sebagai single-hormone AP: hanya insulin, tanpa glukagon.

### 11.2 Plant yang dikendalikan

Controller tidak mengganti persamaan T1D. Keluaran pompa AP menjadi input
`u(t)` pada depot insulin Hovorka:

```text
dS1/dt = uAP(t) - S1/tmaxI
dS2/dt = S1/tmaxI - S2/tmaxI
```

Semua persamaan `Q1`, `Q2`, `I`, `x1`, `x2`, dan `x3`, termasuk clipping EGP,
tetap mengikuti Bagian 10.

### 11.3 Model sensor CGM

Controller membaca `Gs`, bukan `G` secara langsung. Model sensor edukatif awal
menggunakan lag orde pertama:

```text
dGs/dt = (GmgdL - Gs)/Ts
```

Keterangan:

- `GmgdL = 18.0182*Q1/VG`;
- `Gs` adalah pembacaan sensor dalam mg/dL; dan
- `Ts` adalah konstanta waktu sensor.

Nilai awal yang akan diuji adalah `Ts = 10 min`, dengan uji sensitivitas pada
`5`, `10`, dan `15 min`. Hovorka dkk. (2004) membahas delay fisiologis sensor
sekitar 10-15 menit. Bentuk lag orde pertama juga digunakan pada implementasi
sensor dalam perangkat lunak GIM Dalla Man/Cobelli.

Noise dan bias sensor sengaja dipisahkan dari persamaan dasar:

```text
SG(n) = Gs(n) + sensorBias + sensorNoise(n)
```

Verifikasi pertama memakai `sensorBias = 0` dan `sensorNoise = 0`. Uji robustness
baru menambahkan keduanya.

### 11.4 PID diskrit

Controller diperbarui setiap `DeltaTc = 1 min`, mengikuti indeks diskrit paper
PID-IFB. Error didefinisikan positif ketika sensor berada di atas target:

```text
e(n) = SG(n) - Gtarget
```

Komponen controller:

```text
P(n) = Kp*e(n)

Ic_raw(n) = Ic(n-1) + (Kp/Ti)*e(n)*DeltaTc

D(n) = Kp*Td*dSGfiltered(n)/dt

PID(n) = P(n) + Ic(n) + D(n)
```

Untuk `DeltaTc = 1 min`, bentuk integral sama dengan persamaan Steil/Ruiz yang
ditulis tanpa faktor waktu eksplisit. Faktor `DeltaTc` tetap dicantumkan agar
implementasi tidak salah jika interval controller berubah.

Derivative memakai sinyal yang difilter, bukan selisih mentah sensor. Filter
awal proyek:

```text
dSGraw(n)      = [SG(n) - SG(n-1)]/DeltaTc
dSGfiltered(n) = alphaD*dSGfiltered(n-1) + (1-alphaD)*dSGraw(n)
```

`alphaD` adalah parameter implementasi yang harus diuji; filter ini bukan
persamaan spesifik dari paper klinis.

### 11.5 Inisialisasi integral dan anti-windup

Saat loop dimulai pada equilibrium:

```text
Ic(0) = uBasal_U_per_hour
P(0)  = 0
D(0)  = 0
```

Dengan demikian, `PID(0)` sudah menghasilkan basal. Basal tidak boleh ditambahkan
lagi di luar PID karena akan terhitung dua kali.

Paper Ruiz dkk. membatasi integral maksimum menjadi tiga kali basal maksimum dan
meresetnya ketika glukosa berada di bawah `60 mg/dL`. Implementasi proyek memakai
dua lapis:

```text
IcMin <= Ic(n) <= IcMax
IcMax = 3*uBasal_U_per_hour
```

Integral dibekukan jika keluaran telah mencapai saturasi dan arah error akan
mendorong keluaran semakin jauh dari rentang pompa. Ini adalah anti-windup
implementasi proyek.

### 11.6 Estimator insulin feedback

PID-IFB memakai estimator insulin tiga-kompartemen dari Ruiz dkk. (2012), dengan
interval diskrit satu menit:

```text
Isc(n) = a11*Isc(n-1) + b1*ID(n-1)

Ip(n) = a21*Isc(n-1) + a22*Ip(n-1) + b2*ID(n-1)

Ieff(n) = a31*Isc(n-1) + a32*Ip(n-1)
          + a33*Ieff(n-1) + b3*ID(n-1)

IFB(n) = gamma1*Isc(n) + gamma2*Ip(n) + gamma3*Ieff(n)
```

Parameter paper:

| Parameter | Nilai |
|---|---:|
| `a11` | 0.9802 |
| `a21` | 0.014043 |
| `a31` | 0.000127 |
| `a22` | 0.98582 |
| `a32` | 0.017889 |
| `a33` | 0.98198 |
| `b1` | 1.1881 |
| `b2` | 0.0084741 |
| `b3` | 0.00005 |
| `gamma1` | 0.64935 |
| `gamma2` | 0.34128 |
| `gamma3` | 0.0093667 |

Keluaran sebelum pembatasan:

```text
gammaSum = gamma1 + gamma2 + gamma3

ID_PID(n) = (1 + gammaSum)*PID(n)

ID_PID_IFB(n) = (1 + gammaSum)*PID(n) - IFB(n)
```

`ID` dan seluruh state estimator harus memakai satuan yang sama dengan parameter
diskrit paper. State estimator diinisialisasi pada equilibrium yang konsisten
dengan basal, bukan nol. IFB ini menggantikan usulan sederhana
`Kifb*max(0, I-Ib)`; usulan sederhana tersebut tidak dipakai pada model final.

Untuk laju basal diskrit `IDb`, equilibrium estimator dihitung berurutan:

```text
Isc_b  = b1*IDb/(1-a11)

Ip_b   = (a21*Isc_b + b2*IDb)/(1-a22)

Ieff_b = (a31*Isc_b + a32*Ip_b + b3*IDb)/(1-a33)
```

Kemudian `IFB_b` dihitung dari ketiga state tersebut. Pemeriksaan bumpless start
harus membuktikan bahwa persamaan keluaran PID-IFB menghasilkan kembali `IDb`.

### 11.7 Batas pompa dan konversi ke plant

Laju pompa final:

```text
IDsafe(n) = clamp(ID_PID_IFB(n), 0, IDmax)
```

Jika `IDsafe` dinyatakan dalam U/jam, input untuk plant Hovorka adalah:

```text
uAP_mU_kg_min = IDsafe*1000/(60*weightKg)
```

`IDmax` belum dianggap final sebelum uji numerik. Controller juga harus mempunyai
rate limit agar perintah pompa tidak berubah secara tidak terbatas dalam satu
interval.

### 11.8 Low-glucose suspend

Versi edukasi menetapkan lapisan keselamatan terpisah:

```text
jika SG <= Gsuspend:
    IDsafe = 0
    integral tidak boleh naik
```

Nilai awal yang akan diuji adalah `Gsuspend = 70 mg/dL`. Angka ini merupakan
keputusan konservatif proyek dan bukan salinan setting controller pada paper
Ruiz, yang mendokumentasikan reset integral di `60 mg/dL`.

Karena plant hanya mempunyai insulin dan tidak mempunyai glukagon atau hormon
kontraregulasi, suspend dapat menghentikan insulin baru tetapi tidak dapat
menghilangkan insulin yang sudah berada pada `S1`, `S2`, atau plasma.

### 11.9 Parameter kandidat dari paper

Ruiz dkk. (2012) memberikan:

```text
Gtarget = 120 mg/dL
Ti      = 150 min
Td      = 75 min pada siang hari
Td      = 40 min pada malam hari
Kp      = IDIR/2250
```

Paper menuliskan `Kp` dalam `U/hour per mg/dL` dan `IDIR` sebagai kebutuhan
insulin harian subjek. Pemetaan `IDIR` dan satuan terhadap basal populasi plant
Hovorka proyek harus diuji secara eksplisit. Nilai `Kp` dari paper adalah titik
awal riset, bukan parameter final yang boleh langsung dipasang tanpa verifikasi.

Nilai `Ti = 150 min` diuji sebagai kandidat literal dari paper, tetapi pada plant
Hovorka proyek menghasilkan hipoglikemia terlambat setelah makanan 50 dan 75 g.
Sweep terbatas `KpScale = 0.25/0.5/0.75/1`, `Ti = 150/300/450 min`, dan
`Td = 40/75 min` memilih konfigurasi proyek:

```text
KpScale = 1
Ti      = 300 min
Td      = 75 min
```

Perubahan `Ti` dari 150 menjadi 300 menit adalah hasil tuning numerik proyek,
bukan nilai yang dilaporkan paper. Peralihan parameter siang/malam ditunda.

### 11.10 Klasifikasi mode simulator

| Mode | Informasi makan ke controller | Insulin |
|---|---|---|
| T1D manual | Tidak ada controller | Basal dan bolus pengguna |
| AP PID | Tidak ada | Otomatis dari PID |
| AP PID-IFB | Tidak ada | Otomatis dari PID dengan insulin feedback |
| Hybrid AP | Karbohidrat diumumkan | Bolus/feedforward makan + PID-IFB |

Tahap pertama membandingkan AP PID dan AP PID-IFB terhadap T1D manual. Hybrid AP
menjadi tahap berikutnya karena memerlukan aturan bolus makan tambahan.

### 11.11 Verifikasi wajib

1. **Equilibrium:** tanpa makanan, `SG = Gtarget` dan keluaran sama dengan basal.
2. **Bumpless start:** mengaktifkan AP pada equilibrium tidak membuat lonjakan
   insulin.
3. **Konversi satuan:** U/jam controller sama dengan mU/kg/menit pada plant.
4. **PID terms:** `P`, `I`, dan `D` diuji secara terpisah dengan input sintetis.
5. **Estimator IFB:** equilibrium dan respons impuls cocok dengan persamaan
   diskrit paper.
6. **Makan tidak diumumkan:** uji 25, 50, dan 75 g tanpa bolus manual.
7. **Perbandingan:** PID-IFB tidak boleh memberi insulin kumulatif lebih besar
   dari PID ketika estimasi insulin aktif tinggi.
8. **Anti-windup:** integral tetap terbatas ketika pompa tersaturasi.
9. **Suspend:** insulin menjadi nol pada ambang rendah dan tidak menjadi negatif.
10. **Sensor:** uji delay 5/10/15 menit, bias, noise, dan dropout.
11. **Robustness plant:** sensitivitas insulin dan absorpsi makanan divariasikan.
12. **Numerik:** hasil stabil terhadap pengecilan langkah integrasi plant dan
    interval controller.
13. **Metrik edukasi:** puncak, nadir, AUC, total insulin, waktu 70-180, di bawah
    70, di bawah 54, dan di atas 180 mg/dL.

Hasil pelaksanaan checklist dicatat pada Bagian 11.13.

### 11.12 Sumber peer-reviewed

1. Steil, G.M., Rebrin, K., Darwin, C., Hariri, F., & Saad, M.F. (2006).
   *Feasibility of automating insulin delivery for the treatment of type 1
   diabetes.* Diabetes, 55(12), 3344-3350.
   DOI: https://doi.org/10.2337/db06-0419
2. Steil, G.M., Palerm, C.C., Kurtz, N., et al. (2011). *The effect of insulin
   feedback on closed loop glucose control.* Journal of Clinical Endocrinology
   & Metabolism, 96(5), 1402-1408.
   DOI: https://doi.org/10.1210/jc.2010-2578
3. Ruiz, J.L., Sherr, J.L., Cengiz, E., et al. (2012). *Effect of insulin
   feedback on closed-loop glucose control: a crossover study.* Journal of
   Diabetes Science and Technology, 6(5), 1123-1130.
   DOI: https://doi.org/10.1177/193229681200600517
4. Ly, T.T., Keenan, D.B., Roy, A., et al. (2016). *Automated overnight
   closed-loop control using a proportional-integral-derivative algorithm with
   insulin feedback in children and adolescents with type 1 diabetes at diabetes
   camp.* Diabetes Technology & Therapeutics, 18(6), 377-384.
   DOI: https://doi.org/10.1089/dia.2015.0431
5. Pinsker, J.E., Lee, J.B., Dassau, E., et al. (2016). *Randomized crossover
   comparison of personalized MPC and PID control algorithms for the artificial
   pancreas.* Diabetes Care, 39(7), 1135-1142.
   DOI: https://doi.org/10.2337/dc15-2344

### 11.13 Hasil verifikasi numerik PID-IFB v1

Verifier dapat dijalankan dengan:

```text
node verification/ap-pid-ifb-v1.js
```

Asumsi utama:

```text
Berat badan       = 70 kg
Target            = 120 mg/dL
Update controller = 1 min
Plant RK4         = dt 0.1 min
Sensor lag        = 10 min
Kp                = 0.00401943 U/hour per mg/dL
Ti                = 300 min (hasil tuning proyek)
Td                = 75 min
Basal equilibrium = 0.37682 U/hour
Suspend           = 70 mg/dL
```

Seluruh 13 pemeriksaan struktur dan empat gerbang performa lulus:

| Kelompok | Hasil |
|---|---|
| Equilibrium dan bumpless start | Deviasi maksimum PID-IFB `1.67e-16` |
| Konversi satuan | Lolos |
| Komponen P/I/D sintetis | Lolos |
| Fixed point estimator IFB | Galat `5.55e-17` |
| Dosis-respons makanan | Lolos untuk 25/50/75 g |
| IFB membatasi insulin aktif | Lolos |
| Anti-windup dan suspend | Lolos |
| Delay sensor 5/10/15 menit | Stabil; tidak ada glukosa di bawah 70 mg/dL |
| Sensitivitas insulin 0.8/1.0/1.2 | Stabil; tidak ada glukosa di bawah 70 mg/dL |
| Konvergensi RK4 | Hasil `dt=0.1` konsisten dengan `dt=0.05` |

Respons fully closed-loop PID-IFB tanpa bolus atau informasi makan ke controller:

| Makanan | Puncak | Nadir | Time-in-range 70-180 | Waktu <70 | Nilai akhir |
|---:|---:|---:|---:|---:|---:|
| 25 g | 195.92 mg/dL | 97.95 mg/dL | 93.92% | 0 min | 121.40 mg/dL |
| 50 g | 262.55 mg/dL | 82.38 mg/dL | 89.70% | 0 min | 122.56 mg/dL |
| 75 g | 326.20 mg/dL | 74.71 mg/dL | 88.70% | 0 min | 122.84 mg/dL |

Perbandingan makanan 50 g:

| Controller | Puncak | Nadir | Waktu <70 | Total insulin simulasi |
|---|---:|---:|---:|---:|
| PID | 259.25 mg/dL | 57.24 mg/dL | 119.1 min | 10.70 U |
| PID-IFB | 262.55 mg/dL | 82.38 mg/dL | 0 min | 10.40 U |

IFB sedikit menaikkan puncak, tetapi menghilangkan hipoglikemia terlambat pada
skenario nominal dengan mengurangi pemberian insulin ketika estimasi insulin
aktif masih tinggi. Ini sesuai dengan tujuan mekanisme IFB.

#### Temuan dan keputusan implementasi

1. Koefisien estimator pada paper dicetak dengan presisi terbatas. Jika integral
   PID-IFB diinisialisasi tepat sama dengan basal, pembulatan tersebut menghasilkan
   bump kecil. Implementasi memakai equilibrium eksak:

   ```text
   Ic0 = (IDb + IFB_b)/(1 + gammaSum)
   ```

   Untuk PID tanpa IFB:

   ```text
   Ic0 = IDb/(1 + gammaSum)
   ```

2. Input estimator adalah jumlah insulin yang diberikan selama satu interval
   satu menit. Jika laju pompa memakai U/jam:

   ```text
   doseToEstimator = deliveredRateUph/60
   ```

   Mengirim U/jam langsung ke koefisien estimator akan membuat IFB salah sekitar
   60 kali.

3. Nilai literal `Ti = 150 min` dari paper lolos pemeriksaan struktur tetapi
   gagal gerbang performa pada plant proyek. Untuk makanan 50 g, nadir sekitar
   `65.29 mg/dL` dengan sekitar `105.4 min` di bawah 70. Karena itu `Ti = 300 min`
   dipakai untuk simulator edukasi.

4. Fully closed-loop masih menghasilkan puncak tinggi untuk makanan besar karena
   controller baru bereaksi setelah sensor melihat kenaikan glukosa dan insulin
   subkutan memiliki delay. Hybrid AP dengan informasi makan dapat diteliti
   kemudian untuk memperbaiki puncak postprandial.

### 11.14 Status implementasi AP

- Model dan verifier PID-IFB v1 tersedia dan lolos seluruh uji numerik yang
  tercantum.
- Parameter `Ti = 300 min` adalah tuning khusus plant proyek.
- Hasil belum merupakan validasi klinis atau bukti keselamatan perangkat medis.
- PID-IFB belum dimasukkan ke `script.js` atau UI.
- Model siap dipakai sebagai dasar implementasi simulator edukasi, dengan hasil
  browser tetap harus dibandingkan terhadap verifier ini.

---

## 12. Rancangan Mode Diabetes Tipe 2 v1

> **Status arsip:** model `G-I-A` pada bagian ini tetap disimpan sebagai hasil
> riset dan pembanding, tetapi tidak lagi menjadi model T2D utama aplikasi.
> Implementasi UI memakai model `G-X-I` berbasis Mode Normal pada Bagian 13.

### 12.1 Tujuan dan keputusan struktur

Mode T2D v1 ditujukan untuk menunjukkan respons makan pada subjek virtual dengan
resistansi insulin, respons sel-beta yang melemah, efek incretin yang berkurang,
dan supresi glukagon yang tidak sekuat kondisi sehat. Rentang waktu model adalah
menit sampai jam. Progresi penyakit selama bertahun-tahun tidak dimodelkan.

Kajian awal di `sumber/sumbertipe2/model-t2d-glucose-insulin.md` mengusulkan model
hybrid empat-state `G-X-I-A`. Struktur tersebut tidak dipakai karena persamaan
glukosanya menggunakan `a1*I` secara langsung sehingga state `X` tidak memiliki
pengaruh. Persamaan insulin dan glukagon pada rancangan itu juga tidak membentuk
equilibrium dengan kondisi awal yang dicantumkan.

Keputusan v1 adalah memakai tiga state yang konsisten dengan struktur minimal
Subramanian dkk. (2024):

```text
G(t) = glukosa plasma, mg/dL
I(t) = insulin plasma, microU/mL
A(t) = glukagon plasma, pM
```

Semua jalur regulasi dipusatkan terhadap nilai puasa. Centering ini adalah
modifikasi proyek untuk menghasilkan equilibrium yang dapat digunakan berulang
dalam simulator kontinu; bentuk centered ini bukan persamaan literal paper.

### 12.2 Persamaan T2D v1

```text
dG/dt = -SG*(G-Gb)
        -a1*(I-Ib)*G
        +a2*(A-Ab)
        +Ra(t)

dI/dt = -n1*(I-Ib)
        +betaG*[psi(G)-psi(Gb)]
        +betaE*E(t)*psi(G)

dA/dt = -n2*(A-Ab)
        +betaA*[phi(G)-phi(Gb)]
```

Fungsi nonlinear:

```text
psi(G) = 1.5*(G/K)^h / [1 + (G/K)^h]

phi(G) = exp[-kSuppression*G/18.0182]
```

Pembagian `G/18.0182` pada `phi` mengubah glukosa dari mg/dL menjadi mmol/L,
karena `kSuppression` bersatuan L/mmol. Pada titik puasa
`(G,I,A)=(Gb,Ib,Ab)`, tanpa makanan, seluruh selisih bernilai nol sehingga:

```text
dG/dt = dI/dt = dA/dt = 0
```

Makna tiap jalur:

| Suku | Makna |
|---|---|
| `SG*(G-Gb)` | Pemakaian glukosa yang tidak bergantung insulin |
| `a1*(I-Ib)*G` | Pemakaian glukosa yang bergantung insulin; `a1` rendah merepresentasikan resistansi |
| `a2*(A-Ab)` | Pengaruh perubahan glukagon terhadap keluaran glukosa hati |
| `betaG*[psi(G)-psi(Gb)]` | Sekresi insulin akibat kenaikan glukosa |
| `betaE*E(t)*psi(G)` | Potensiasi sekresi insulin oleh surrogate incretin setelah makan |
| `betaA*[phi(G)-phi(Gb)]` | Perubahan sekresi glukagon akibat glukosa |

### 12.3 Input makanan dan surrogate incretin

Kemunculan glukosa memakai kernel dua-kompartemen Hovorka yang sama dengan mode
T1D, lalu dikonversi menjadi mg/dL/menit:

```text
DG = carbsGram*1000/(180.16*weightKg)                 mmol/kg

UG(t) = DG*AG*age*exp(-age/tmaxG)/(tmaxG^2)          mmol/kg/min

Ra(t) = UG(t)/VG*18.0182                             mg/dL/min
```

Efek incretin v1 memakai sinyal empiris yang dipicu oleh makanan:

```text
E(t) = (carbsGram/50)*(age/tauE)*exp(1-age/tauE)
```

Sinyal tersebut tidak memiliki satuan dan mencapai nilai satu pada `age=tauE`
untuk makanan 50 g. Jika ada beberapa makanan, kontribusinya dijumlahkan.

Paper Subramanian menggunakan profil GLP-1 atau GIP terukur sebagai fungsi input.
Karena simulator tidak mempunyai data hormon terukur, `E(t)` adalah surrogate
proyek, bukan salinan literal model incretin paper. Parameter `betaE` dan `tauE`
harus diperlakukan sebagai parameter kalibrasi edukatif.

### 12.4 Parameter v1

```text
weightKg       = 70 kg
Gb             = 117 mg/dL
Ib             = 12 microU/mL
Ab             = 20 pM
SG             = 0.014 1/min
a1             = 6.6e-5 (microU/mL*min)^-1
a2             = 0.16 mg/dL/(pM*min)
n1             = 0.14 1/min
n2             = 0.08 1/min
betaG          = 3.2 microU/mL/min
betaE          = 0.55 microU/mL/min
betaA          = 5.0 pM/min
K              = 306 mg/dL
h              = 1.27
kSuppression   = 0.15 L/mmol
AG             = 0.8
VG             = 0.16 L/kg
tmaxG          = 40 min
tauE           = 35 min
```

Asal dan status parameter:

| Parameter | Status |
|---|---|
| `SG`, `n1`, `n2`, `K`, `h` | Nilai populasi/fixed dari Subramanian dkk. (2024) |
| `a1` | Rata-rata T2D Subramanian yang dikonversi dari insulin `10 pM` ke `microU/mL` dengan `1 microU/mL ~= 6 pM` |
| `betaA`, `kSuppression` | Nilai rata-rata T2D Subramanian; dipakai dalam formulasi centered proyek |
| `a2=0.16` | Koefisien glukagon efektif hasil kalibrasi proyek; nilai rata-rata paper `0.26` terlalu dominan setelah formulasi diubah menjadi centered |
| `AG`, `VG`, `tmaxG`, `weightKg` | Parameter input makanan dari plant Hovorka proyek |
| `Gb=117` | Profil subjek virtual T2D awal/ringan; konsisten dengan nilai jaringan yang dikutip Banzi dkk. dari Vahidi |
| `Ib=12`, `Ab=20` | Pilihan kondisi basal subjek virtual; bukan rata-rata pasien universal |
| `betaG`, `betaE`, `tauE` | Kalibrasi proyek untuk respons makan v1 |

Nilai `117 mg/dL` tidak dipilih sebagai batas diagnosis. Ia merupakan baseline
subjek virtual edukatif. Mode ini tidak boleh digunakan untuk mengklasifikasikan
atau mendiagnosis pengguna.

### 12.5 Sumber peer-reviewed

1. Subramanian, V., Bagger, J.I., Harihar, V., Holst, J.J., Knop, F.K., &
   Vilsboll, T. (2024). *An extended minimal model of OGTT: estimation of alpha-
   and beta-cell dysfunction, insulin resistance, and the incretin effect.*
   American Journal of Physiology-Endocrinology and Metabolism, 326, E182-E205.
   DOI: https://doi.org/10.1152/ajpendo.00278.2023
2. Lopez-Palau, N.E., & Olais-Govea, J.M. (2020). *Mathematical model of blood
   glucose dynamics by emulating the pathophysiology of glucose metabolism in
   type 2 diabetes mellitus.* Scientific Reports, 10, 12697.
   DOI: https://doi.org/10.1038/s41598-020-69629-0
3. Banzi, W., Kambutse, I., Dusabejambo, V., et al. (2021). *Mathematical
   Modelling of Glucose-Insulin System and Test of Abnormalities of Type 2
   Diabetic Patients.* International Journal of Mathematics and Mathematical
   Sciences, 2021, 6660177. DOI: https://doi.org/10.1155/2021/6660177
4. Yang, B., Li, J., Haller, M.J., Schatz, D.A., & Rong, L. (2023). *Modeling
   the progression of Type 2 diabetes with underlying obesity.* PLOS
   Computational Biology, 19(2), e1010914.
   DOI: https://doi.org/10.1371/journal.pcbi.1010914

Subramanian adalah sumber struktur utama untuk dinamika episode OGTT. Lopez-Palau
dan Banzi mendukung jalur patofisiologi dan perilaku T2D. Yang dkk. dipakai hanya
untuk konteks progresi jangka panjang dan tidak menjadi sumber parameter episode
makan v1.

### 12.6 Verifikasi numerik v1

Verifier dapat dijalankan dengan:

```text
node verification/t2d-v1.js
```

Integrator menggunakan RK4. Pemeriksaan yang sudah lulus:

1. equilibrium puasa 1000 menit tanpa drift;
2. pemulihan ke equilibrium setelah perturbasi terpisah pada `G`, `I`, dan `A`;
3. positivitas seluruh state untuk makanan sampai 125 g;
4. dosis-respons glukosa monoton untuk 25, 50, dan 75 g;
5. insulin meningkat dan glukagon tertekan setelah makanan;
6. peningkatan sensitivitas insulin menurunkan puncak glukosa;
7. respons tidak bergantung pada waktu absolut makanan;
8. konvergensi RK4 pada `dt=1`, `0.5`, `0.25`, dan `0.1 min`.

Verifikasi tambahan tanpa clipping state juga lulus. Seluruh state tetap positif
pada skenario nominal dan stress 125 g, sehingga positivitas bukan hasil pemaksaan
`max(0,state)` oleh integrator.

Neraca input makanan 50 g:

```text
Integral analitik Ra = 357.18646917 mg/dL
Integral numerik Ra  = 357.18646997 mg/dL
Galat                 = 8.01e-7 mg/dL
```

Linearisasi numerik di equilibrium menghasilkan matriks Jacobian:

```text
[-0.014000  -0.007722   0.160000]
[ 0.009164  -0.140000   0.000000]
[-0.015716   0.000000  -0.080000]
```

Eigenvalue independen yang dihitung untuk audit adalah kira-kira
`-0.04721 +/- 0.03817i` dan `-0.13958`. Semua bagian real negatif. Verifier
juga memeriksa syarat Routh-Hurwitz kubik secara langsung dan menyatakan
equilibrium stabil lokal.

Hasil nominal `dt=0.1 min`, makanan pada menit ke-60:

| Karbohidrat | Puncak glukosa | Puncak insulin | Nadir glukagon | Glukosa akhir 1000 min |
|---:|---:|---:|---:|---:|
| 25 g | 153.86 mg/dL | 14.97 microU/mL | 14.02 pM | 117.0000 mg/dL |
| 50 g | 194.32 mg/dL | 18.17 microU/mL | 9.16 pM | 117.0000 mg/dL |
| 75 g | 238.28 mg/dL | 21.50 microU/mL | 5.38 pM | 117.0000 mg/dL |

Untuk makanan 50 g, puncak glukosa pada `dt=1`, `0.5`, `0.25`, dan `0.1 min`
berturut-turut adalah `194.31771`, `194.31771`, `194.31774`, dan
`194.31801 mg/dL`. Selisih maksimum kurang dari `0.001 mg/dL`.

Sweep `a1` sebesar `0.5x`, `1x`, dan `2x` menghasilkan puncak berturut-turut
`195.19`, `194.32`, dan `192.63 mg/dL`. Arah respons sudah benar.

Audit kontribusi jalur untuk makanan 50 g menemukan:

| Perubahan | Puncak glukosa | Selisih dari nominal |
|---|---:|---:|
| Nominal T2D | 194.32 mg/dL | - |
| Aksi insulin dihilangkan (`a1=0`) | 196.09 mg/dL | +1.77 mg/dL |
| `a1` kontrol paper (`3.0e-4`) | 188.66 mg/dL | -5.66 mg/dL |
| Kontribusi glukagon ke glukosa dihilangkan (`a2=0`) | 248.90 mg/dL | +54.58 mg/dL |
| Surrogate incretin dihilangkan | 194.90 mg/dL | +0.59 mg/dL |

Nilai literal rata-rata paper `a2=0.26` membuat jalur glukagon terlalu dominan
dalam formulasi centered: perbedaan sensitivitas insulin T2D-versus-kontrol hanya
`3.56 mg/dL`. Sweep `a2=0.05-0.26` memilih `a2=0.16` sebagai kompromi pertama:
pemisahan sensitivitas insulin menjadi `5.66 mg/dL`, respons 50 g tetap berada
pada rentang hiperglikemia edukatif, dan stress 125 g tetap positif tanpa clipping.
Verifier memakai gerbang proyek minimal `5 mg/dL` untuk memastikan resistansi
insulin terlihat secara material. Gerbang ini sekarang **lulus**. Ambang tersebut
adalah kriteria implementasi proyek, bukan ambang klinis atau nilai dari paper.

### 12.7 Status dan batas penggunaan

- Struktur matematis v1, equilibrium, stabilitas lokal, neraca makanan,
  positivitas tanpa clipping, dan pemeriksaan numerik dasar sudah lulus.
- Gerbang atribusi fisiologis lulus setelah koefisien glukagon efektif `a2`
  dikalibrasi dari nilai literal paper `0.26` menjadi nilai proyek `0.16`;
  `readyForUI = true` pada verifier.
- Model belum difit ke time-series satu pasien atau rata-rata kohort OGTT.
- Surrogate incretin belum divalidasi terhadap data GLP-1/GIP.
- Nilai `betaG`, `betaE`, dan `tauE` adalah hasil kalibrasi proyek, bukan nilai
  langsung dari paper.
- Hysteresis glukagon `k1/k2` dari Subramanian belum dimasukkan pada v1.
- Model belum diintegrasikan ke `script.js` atau UI.
- Status yang tepat adalah **lulus struktural dan lulus gerbang fungsional proyek**.
  Model dapat menjadi base implementasi UI edukatif, tetapi belum tervalidasi
  terhadap time-series OGTT dan bukan model T2D tervalidasi secara klinis.

---

## 13. Model T2D Utama Berbasis Mode Normal

### 13.1 Tujuan

Model T2D utama memakai persamaan `G-X-I` Mode Normal tanpa menambah state atau
controller baru. Tujuannya adalah menyediakan perbandingan Normal-versus-T2D
yang mudah dipahami karena kedua mode menerima input makan dan memakai struktur
matematis yang sama.

Persamaannya tetap:

```text
dG/dt = -p1*(G-Gb) - X*G + D(t)
dX/dt = -p2*X + p3*(I-Ib)
dI/dt = p6*[G-p5]+*tau - n*(I-Ib)
```

PID tidak digunakan. Persamaan sekresi insulin tetap menjadi umpan balik
fisiologis internal model.

### 13.2 Profil parameter T2D

| Parameter | Normal | T2D sederhana | Makna perubahan |
|---|---:|---:|---|
| `Gb` | 92 mg/dL | 117 mg/dL | Baseline glukosa subjek virtual lebih tinggi |
| `Ib` | 7.3 microU/mL | 12 microU/mL | Insulin puasa kompensatorik lebih tinggi |
| `p1` | 0.03082 | 0.024 1/min | Glucose effectiveness lebih rendah |
| `p2` | 0.02093 | 0.02093 1/min | Dipertahankan untuk isolasi perubahan utama |
| `p3` | 1.062e-5 | 2.655e-6 | Sensitivitas insulin 25% dari Normal |
| `n` | 0.3 | 0.3 1/min | Clearance insulin dipertahankan |
| `p6` | 0.003349 | 0.001675 | Respons sel-beta 50% dari Normal |
| `p5` | 94 mg/dL | 119 mg/dL | Ambang sekresi tetap di atas `Gb` agar equilibrium valid |

Persentase `p3` dan `p6`, serta nilai `p1`, adalah **pengaturan tingkat keparahan
proyek**. Nilai tersebut belum merupakan hasil fitting kohort T2D dan tidak boleh
ditafsirkan sebagai nilai pasien.

### 13.3 Verifikasi numerik

Verifier:

```text
node verification/t2d-normal-derived-v1.js
```

Hasil dengan input makan Mode Normal yang sama:

| Input | Normal: puncak G | T2D sederhana: puncak G | Kenaikan dari baseline T2D | Puncak insulin T2D |
|---|---:|---:|---:|---:|
| Kecil (`A=80`) | 108.77 mg/dL | 138.55 mg/dL | 21.55 mg/dL | 20.95 microU/mL |
| Sedang (`A=150`) | 122.54 mg/dL | 157.05 mg/dL | 40.05 mg/dL | 29.80 microU/mL |
| Besar (`A=220`) | 135.94 mg/dL | 175.37 mg/dL | 58.37 mg/dL | 38.48 microU/mL |

Seluruh pemeriksaan berikut lulus:

1. equilibrium Normal dan T2D tidak drift selama 1000 menit;
2. perturbasi terpisah pada `G`, `X`, dan `I` kembali ke equilibrium;
3. seluruh state tetap nonnegatif;
4. respons makanan kecil, sedang, dan besar bersifat monoton;
5. kenaikan glukosa T2D lebih besar daripada Normal untuk input identik;
6. mengembalikan `p3` atau `p6` ke nilai Normal menurunkan puncak glukosa;
7. respons tidak bergantung pada waktu absolut makanan;
8. hasil konvergen pada `dt=1`, `0.5`, `0.25`, dan `0.1 min`.

Untuk input sedang, puncak pada `dt=1` dan `dt=0.1` masing-masing adalah
`157.06794` dan `157.05463 mg/dL`, dengan selisih sekitar `0.013 mg/dL`.

### 13.4 Keputusan implementasi

| Model | Kegunaan utama | Batas utama |
|---|---|---|
| T2D utama `G-X-I` | Perbandingan langsung dengan Normal, plant untuk eksperimen kontrol, dan penjelasan resistansi insulin | Tidak mempunyai state glukagon atau incretin |
| T2D riset `G-I-A` Bagian 12 | Eksplorasi insulin, glukagon, dan surrogate incretin | Diarsipkan sebagai pembanding dan tidak dipakai UI |

Model `G-X-I` ditetapkan sebagai **plant T2D utama** dan sudah dipasang pada UI.
Struktur persamaan sama dengan Mode Normal, sedangkan `Gb`, `Ib`, `p1`, `p3`,
`p5`, dan `p6` memakai profil T2D. Pemisahan ini memungkinkan controller P/PI/PID
ditambahkan kemudian sebagai input eksternal tanpa mengubah parameter plant.

Model `G-I-A` Bagian 12 tidak dihapus agar keputusan dan hasil riset sebelumnya
tetap dapat diaudit. Model tersebut tidak lagi dijalankan oleh `script.js`.

### 13.5 Posisi controller T2D

Plant T2D tetap memakai parameter pada Bagian 13.2. Controller nantinya menjadi
lapisan eksternal dengan kandidat setpoint awal `100 mg/dL`; controller tidak
mengubah `Gb`, `p1`, `p3`, atau `p6`. Tanpa controller, plant kembali ke baseline
T2D `117 mg/dL`. Karena itu UI menyebut kondisi basal tersebut **di atas target**,
bukan kondisi Normal. Struktur aktuator insulin, gain, saturasi, anti-windup, dan
uji keselamatan perlu dirancang serta diverifikasi sebelum mode kontrol T2D
diaktifkan.

---

## 14. Rancangan PI untuk Plant T2D G-X-I

### 14.1 Tujuan dan arsitektur

Controller PI v1 digunakan untuk eksperimen sistem kontrol: mempertahankan plant
T2D pada kandidat target `100 mg/dL` tanpa mengganti parameter T2D Bagian 13.2.

```text
target 100 mg/dL -> error -> PI -> saturasi/suspend -> insulin subkutan
       ^                                                   |
       |                                                   v
       +---- sensor glukosa <- plant T2D G-X-I <- S1/S2 ---+
```

Ini merupakan algoritma proyek untuk simulator pendidikan. Ia bukan algoritma
pompa insulin tervalidasi dan keluarannya tidak boleh dipakai sebagai rekomendasi
dosis.

### 14.2 Aktuator insulin eksternal

Plant T2D diperluas dengan dua depot insulin subkutan dan satu sensor lag:

```text
dS1/dt = u - S1/tmaxI
dS2/dt = S1/tmaxI - S2/tmaxI

dI/dt = p6*[G-p5]+*tau - n*(I-Ib)
        + S2/(tmaxI*VI)

dGs/dt = (G-Gs)/tauSensor
```

`u` dikonversi dari keluaran controller U/jam menjadi mU/kg/menit:

```text
u = rateUph*1000/(60*weightKg)
```

Struktur dua depot dan nilai kandidat `VI=0.12 L/kg`, `tmaxI=55 min` mengikuti
kompartemen absorpsi insulin Hovorka yang sudah digunakan pada plant T1D proyek.
Penggabungannya dengan persamaan sekresi endogen `G-X-I` adalah adaptasi proyek.

### 14.3 Persamaan PI

Konvensi error dibuat positif ketika glukosa berada di atas target:

```text
e(k) = Gs(k) - Gtarget

rateRaw(k) = ubias + Kp*e(k) + Ki*J(k)
Ki         = Kp/Ti

rate(k) = clamp(rateRaw, 0, umax)
```

Integral diskrit:

```text
J(k+1) = J(k) + e(k)*Ts
```

Integral dibekukan jika keluaran tersaturasi dan error mendorong lebih jauh ke
arah saturasi. Laju insulin menjadi nol jika glukosa sensor berada pada atau di
bawah ambang suspend.

Parameter kandidat hasil tuning numerik proyek:

```text
Gtarget   = 100 mg/dL
Ts        = 1 min
tauSensor = 10 min
Kp        = 0.035 U/hour per mg/dL
Ti        = 300 min
Ki        = 0.000116667 U/hour per (mg/dL*min)
umax      = 10 U/hour
Gsuspend  = 75 mg/dL
```

Nilai gain tersebut bukan hasil paper atau parameter terapi pasien. Paper PI/PID
pada Bagian 11 mendukung konsep feedback insulin, tetapi parameter T2D v1 dituning
khusus untuk satu plant virtual proyek.

### 14.4 Bias dan equilibrium target

Pada target tanpa makanan, sekresi endogen bernilai nol karena `Gtarget < p5`.
Equilibrium dihitung dari persamaan plant, bukan ditebak:

```text
Xtarget = p1*(Gb-Gtarget)/Gtarget
Itarget = Ib + p2*Xtarget/p3
ubias   = n*(Itarget-Ib)*VI
```

Untuk profil T2D v1 dan berat 70 kg:

```text
Xtarget = 0.00408 1/min
Itarget = 44.1636 microU/mL
ubias   = 1.15789 mU/kg/min
        = 4.86314 U/hour
```

Nilai bias tersebut tinggi karena plant mempunyai suku `-p1*(G-Gb)` yang selalu
mendorong glukosa kembali ke baseline `117 mg/dL`. Ini adalah konsekuensi model
minimal ketika dipaksa menetap pada `100 mg/dL`, bukan rekomendasi basal nyata.

### 14.5 Verifikasi numerik PI v1

Verifier:

```text
node verification/t2d-pi-v1.js
```

Hasil utama `dt=0.1 min`:

| Skenario | Puncak G | Nadir G | G akhir | Laju maksimum | Waktu <70 |
|---|---:|---:|---:|---:|---:|
| Tanpa kontrol, tanpa makan | 117.00 | 117.00 | 117.00 | 0 | 0 min |
| PI, tanpa makan | 117.00 | 98.95 | 99.30 | 5.58 U/jam | 0 min |
| PI + makan kecil | 118.82 | 97.98 | 99.02 | 5.97 U/jam | 0 min |
| PI + makan sedang | 135.98 | 97.08 | 98.64 | 6.65 U/jam | 0 min |
| PI + makan besar | 152.87 | 96.26 | 98.29 | 7.32 U/jam | 0 min |

Pemeriksaan berikut lulus:

1. equilibrium terkontrol di `100 mg/dL` tidak drift;
2. dari baseline T2D `117 mg/dL`, PI mendekati target tanpa turun di bawah 70;
3. dosis-respons makanan kecil/sedang/besar monoton;
4. seluruh state tetap nonnegatif;
5. saturasi atas `10 U/jam` bekerja;
6. conditional integration menghentikan windup pada saturasi;
7. suspend pada `75 mg/dL` menghasilkan keluaran nol;
8. hasil konsisten pada `dt=0.5`, `0.25`, dan `0.1 min`.

### 14.6 Status implementasi

- Plant dan verifier PI v1 tersedia dan seluruh gerbang numerik di atas lulus.
- Tab `T2D + PI` sudah menggunakan plant G-X-I, aktuator dua depot, sensor lag,
  controller PI, serta input makanan yang sama dengan Mode Normal/T2D.
- Regression test aplikasi membandingkan respons nominal PI terhadap verifier:
  puncak makanan sedang `135.98 mg/dL` dan glukosa akhir `98.64 mg/dL`.
- Slider resistansi `0.5-2.0` mengubah `p3` efektif; bias equilibrium controller
  dihitung ulang untuk nilai yang dipilih. Skenario tanpa makan dan makan sedang
  pada kedua batas diuji sampai 2400 menit, kembali dalam `2 mg/dL` dari target
  tanpa glukosa di bawah `70 mg/dL`. Perubahan slider dari `1x` ke `2x` pada
  menit ke-500 juga diuji; glukosa akhir `99.61 mg/dL`, nadir `98.95 mg/dL`.
  Perubahan `2x` ke `0.5x` menghasilkan nadir `79.92 mg/dL` dan pulih ke
  `100.70 mg/dL` pada menit ke-2400. Insulin di depot membuat respons awal
  sesudah perubahan parameter tetap berlangsung walaupun bias dihitung ulang.
- Bias `4.86 U/jam` dan insulin target `44.16 microU/mL` menunjukkan bahwa angka
  keluaran tidak layak ditampilkan sebagai saran dosis dunia nyata.
- UI memberi label **sinyal kontrol model** dan peringatan pendidikan.
- Model dapat digunakan untuk eksperimen tugas kontrol dan perbandingan tanpa
  kontrol versus PI, tetapi belum mempunyai validasi klinis atau robustness
  antar-populasi.

---

*Dokumen ini memuat Mode Normal, T1D basal-bolus, Artificial Pancreas PID-IFB,
arsip riset T2D G-I-A, dan model T2D utama berbasis Mode Normal. Seluruh mode
hanya untuk pendidikan dan bukan prediksi, diagnosis, atau rekomendasi klinis.*
