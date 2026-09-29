# Model Glukosa-Insulin: Mode Diabetes Tipe 2 (T2D)

> **ARSIP RISET AWAL:** persamaan hybrid `G-X-I-A` dalam dokumen ini tidak dipakai
> sebagai base implementasi karena state `X` tidak terhubung ke persamaan glukosa,
> satuan belum konsisten, dan kondisi awalnya bukan equilibrium. Spesifikasi T2D
> v1 yang telah melewati verifikasi numerik dasar berada pada Bagian 12 dokumen
> `../model-normal-glucose-insulin.md`. Dokumen ini dipertahankan sebagai catatan
> penelusuran paper dan rancangan yang ditinggalkan.

Dokumentasi ini merangkum model matematis, sumber, parameter, modifikasi, dan rancangan
verifikasi untuk **Mode Diabetes Tipe 2 (T2D)** pada simulator sistem glukosa-insulin.
Dokumen ini dibangun di atas fondasi Mode Normal (Bergman–Pacini–Friis-Jensen) dan
dirancang agar kompatibel dengan simulator yang sama.

> **Status:** spesifikasi riset; belum melewati verifikasi numerik penuh. Siap untuk
> implementasi bertahap dengan verifikasi wajib di setiap lapisan.

---

## 1. Ringkasan Model

Mode T2D merepresentasikan kondisi tubuh penderita diabetes tipe 2, yang ditandai oleh
tiga abnormalitas fisiologis utama:

1. **Resistansi insulin** — sensitivitas jaringan perifer terhadap insulin menurun,
   sehingga insulin tidak efektif mendorong pengambilan glukosa.
2. **Disfungsi sel-β (β-cell dysfunction)** — pankreas masih mensekresi insulin,
   tetapi kapasitasnya berkurang, terutama respons fase pertama dan efek incretin.
3. **Disfungsi sel-α (α-cell dysfunction)** — supresi sekresi glukagon oleh glukosa
   terganggu, sehingga glukagon tetap tinggi meski glukosa darah sudah naik.

Berbeda dengan Mode Normal yang merupakan **closed-loop fisiologis sempurna** (pankreas
merespons glukosa secara alami), Mode T2D mempertahankan loop tersebut tetapi dengan
**parameter yang dimodifikasi** untuk mencerminkan ketiga abnormalitas di atas.

Mode T2D juga menambahkan **glukagon A(t)** sebagai state keempat, karena peran glukagon
dalam patofisiologi T2D tidak dapat diabaikan.

### Pemetaan ke Block Diagram Kontrol

| Blok | Mode Normal | Mode T2D |
|---|---|---|
| Input / Setpoint | `p5` (ambang glukosa) | `p5` (sama) |
| Controller | `dI/dt` — sekresi pankreas normal | `dI/dt` — sekresi melemah + incretin turun |
| Process (Plant) | `dG/dt`, `dX/dt` | `dG/dt` (+ kontribusi glukagon), `dX/dt` (sensitivitas turun) |
| Feedback path | `G(t)` → `[G-p5]⁺` | `G(t)` → `[G-p5]⁺` (tetap, tapi gain berkurang) |
| State tambahan | — | `dA/dt` (dinamika glukagon) |
| Disturbance | `D(t)` (makan) | `D(t)` (sama) |

---

## 2. Persamaan Matematis

Sistem **4 persamaan diferensial nonlinear, closed-loop**:

```
dG(t)/dt = -(SG + a1·I(t))·G(t) + a2·A(t) + D(t)          ... (1)

dX(t)/dt = -p2·X(t) + p3·(I(t) - Ib)                        ... (2)

dI(t)/dt = -n1·(I(t) - Ib) + γ1·ψ(G(t)) + γ3·Φ(G(t))       ... (3)

dA(t)/dt = -n2·(A(t) - Ab) + γ2·φ(G(t) - τ1)                ... (4)
```

**Fungsi nonlinear tambahan:**

```
ψ(G) = 1.5·(G/K)^h / (1 + (G/K)^h)     ← fungsi Hill, sekresi insulin dari glukosa
                                            K = 17.0 mM, h = 1.27 (dari data islet)

φ(G) = exp(-k1·G)                         ← respons glukagon: turun eksponensial saat G naik

Φ(G) = γ3·G(t)                            ← incretin-potentiated insulin secretion
                                            (disederhanakan sebagai fungsi linier terhadap G)
```

**Keterangan variabel:**
- `G(t)` — konsentrasi glukosa plasma (mg/dL)
- `X(t)` — efek insulin aktif di jaringan remote (1/menit); dipertahankan dari Mode Normal
- `I(t)` — konsentrasi insulin plasma (μU/mL)
- `A(t)` — konsentrasi glukagon plasma (pM)
- `D(t)` — input gangguan makan; fungsi identik dengan Mode Normal: `A·t·e^(-t/τ)/τ²`
- `SG` — glucose effectiveness (konstanta, identik dengan p1 di Mode Normal = 0.03082/min)
- `a1` — laju klirens glukosa bergantung insulin (**parameter T2D utama: resistansi insulin**)
- `a2` — efek glukagon terhadap produksi glukosa hepatik
- `γ1` — laju sekresi insulin bergantung glukosa
- `γ2` — amplitudo sekresi glukagon
- `γ3` — laju sekresi insulin potentiated oleh incretin (**parameter T2D: melemah di T2D**)
- `k1` — konstanta supresi glukagon oleh glukosa (**parameter T2D: melemah di T2D**)
- `n1` — klirens insulin plasma
- `n2` — klirens glukagon plasma
- `Ab` — glukagon basal (pM)

**Kondisi awal (baseline/istirahat):**
```
G(0) = Gb = 126 mg/dL    ← fasting hyperglycemia khas T2D (> normal 92 mg/dL)
X(0) = 0
I(0) = Ib                ← kompensasi awal: bisa lebih tinggi dari normal
A(0) = Ab                ← glukagon basal lebih tinggi dari normal
```

---

## 3. Sumber (Rantai Sitasi)

| Bagian Model | Sumber |
|---|---|
| Persamaan glukosa `dG/dt` — adaptasi dari minimal model Bergman | Bergman, R.N., Ider, Y.Z., Bowden, C.R., & Cobelli, C. (1979). Quantitative estimation of insulin sensitivity. *American Journal of Physiology*, 236(6), E667–E677. |
| Persamaan efek insulin `dX/dt` | Bergman et al. (1979) — identik dengan Mode Normal. |
| Persamaan insulin `dI/dt` dengan fungsi Hill dan suku incretin | Subramanian, V., Bagger, J.I., Harihar, V., Holst, J.J., Knop, F.K., & Villsbøll, T. (2024). An extended minimal model of OGTT: estimation of α- and β-cell dysfunction, insulin resistance, and the incretin effect. *American Journal of Physiology — Endocrinology and Metabolism*, 326, E182–E205. DOI: 10.1152/ajpendo.00278.2023 |
| Persamaan glukagon `dA/dt` dan fungsi respons glukagon hiperbolik-eksponensial | Subramanian et al. (2024) — Persamaan (3) dan fungsi φ(G). |
| Fungsi Hill `ψ(G)` untuk insulin dose-response; K = 17.0 mM, h = 1.27 dari data islet | Subramanian et al. (2024) — Persamaan (7); parameter dari data islet Rorsman et al. |
| Kuantifikasi perubahan parameter T2D vs kontrol sehat (`a1`, `k1`, `γ2`, `γ3`) | Subramanian et al. (2024) — Tabel 1, 2, dan 3; N = 8 pasien T2D dan 8 kontrol. |
| Konsep resistansi insulin sebagai penurunan `a1`; homeostasis glukosa berbasis Bergman | Yang, B., Li, J., Haller, M.J., Schatz, D.A., & Rong, L. (2023). Modeling the progression of Type 2 diabetes with underlying obesity. *PLOS Computational Biology*, 19(2), e1010914. DOI: 10.1371/journal.pcbi.1010914 |
| Tiga abnormalitas T2D: resistansi insulin perifer, produksi glukosa hepatik berlebih, disfungsi sel-β | López-Palau, N.E., & Olais-Govea, J.M. (2020). Mathematical model of blood glucose dynamics by emulating the pathophysiology of glucose metabolism in type 2 diabetes mellitus. *Scientific Reports*, 10, 12697. DOI: 10.1038/s41598-020-69629-0 |
| Struktur subsystem 4-kompartemen (jantung-paru, hati, jaringan, pankreas); validasi dengan data pasien T2D Rwanda | Banzi, W., Kambutse, I., Dusabejambo, V., Rutaganda, E., Minani, F., Niyobuhungiro, J., Mpinganzima, L., & Ntaganda, J.M. (2021). Mathematical Modelling of Glucose-Insulin System and Test of Abnormalities of Type 2 Diabetic Patients. *International Journal of Mathematics and Mathematical Sciences*, 2021, 6660177. DOI: 10.1155/2021/6660177 |
| Fungsi input makan `D(t)` | Hovorka, R. et al. (2004). Nonlinear model predictive control of glucose concentration in subjects with type 1 diabetes. *Physiological Measurement*, 25(4), 905–920. — identik dengan Mode Normal. |

**Catatan penting tentang pendekatan hybrid ini:**

Model T2D ini adalah **adaptasi edukatif** yang mengombinasikan elemen dari beberapa paper.
Persamaan (1)–(2) diadaptasi dari Bergman (1979) dengan penambahan suku glukagon `a2·A(t)`.
Persamaan (3) diadaptasi dari Subramanian et al. (2024) yang telah memvalidasinya terhadap
data OGTT pasien T2D. Persamaan (4) merupakan penyederhanaan dari model glukagon
Subramanian et al. (2024). Gabungan keempat persamaan ini **belum divalidasi sebagai
satu sistem tunggal** dalam satu paper, sehingga harus disebut model edukasi hasil adaptasi
dan wajib diverifikasi sebagai sistem baru.

---

## 4. Parameter

### 4.1 Parameter Dipertahankan dari Mode Normal

Parameter berikut **tidak berubah** karena mencerminkan fisiologi dasar yang tidak spesifik
untuk T2D:

| Simbol | Nilai | Satuan | Arti Fisiologis |
|---|---|---|---|
| SG (≡ p1) | 0.03082 | 1/min | Glucose effectiveness |
| p2 | 0.02093 | 1/min | Peluruhan efek insulin aktif |
| p3 | 1.062×10⁻⁵ | L/(min²·mU) | Laju aktivasi insulin di jaringan remote |
| Gb | 126 | mg/dL | Glukosa basal T2D (fasting hyperglycemia) |
| mealTau | 40 | menit | Konstanta waktu absorpsi makanan |
| K (Hill) | 17.0 | mM | Konstanta setengah-maksimum fungsi Hill |
| h (Hill) | 1.27 | — | Koefisien Hill |

**Catatan:** `Gb` dinaikkan dari 92 (normal) menjadi 126 mg/dL sebagai representasi
*fasting hyperglycemia* khas T2D (kriteria diagnostik ADA: FPG ≥ 126 mg/dL).

### 4.2 Parameter T2D Utama (Dimodifikasi dari Kontrol Sehat)

Nilai berikut diturunkan dari Subramanian et al. (2024), Tabel 1–3, perbandingan
7 kontrol sehat vs 8 pasien T2D (75 g OGTT):

| Simbol | Nilai Kontrol Sehat | Nilai T2D | Rasio | Sumber |
|---|---|---|---|---|
| a1 (×10⁻⁴) | 5.0 (10 pM·min)⁻¹ | **1.1** (10 pM·min)⁻¹ | ↓ 4.5× | Subramanian et al. (2024), Tabel 3 |
| k1 (glucagon suppression) | 0.25 mM⁻¹ | **0.15** mM⁻¹ | ↓ 40% | Subramanian et al. (2024), Tabel 3 |
| γ2 (glucagon magnitude) | 2.6 pM·min⁻¹ | **5.0** pM·min⁻¹ | ↑ 1.9× | Subramanian et al. (2024), Tabel 3 |
| γ3 (incretin-dep. secretion) | 0.0068 (mg/dL·min)⁻¹ | **0.0017** (mg/dL·min)⁻¹ | ↓ 4× | Subramanian et al. (2024), Tabel 3 |
| n1 (insulin clearance) | 0.14 min⁻¹ | 0.14 min⁻¹ | = | Subramanian et al. (2024), fixed |
| n2 (glucagon clearance) | 0.08 min⁻¹ | 0.08 min⁻¹ | = | Subramanian et al. (2024), fixed |

### 4.3 Parameter Baru (Glukagon)

| Simbol | Nilai | Satuan | Arti | Sumber |
|---|---|---|---|---|
| Ab | 20.0 | pM | Glukagon basal T2D | Subramanian et al. (2024), Tabel 2 |
| a2 | 0.26 | mg/dL·(pM·min)⁻¹ | Efek glukagon pada HGP | Subramanian et al. (2024), Tabel 2 |
| k2 (recovery) | 0.56 | mM⁻¹ | Konstanta pemulihan glukagon | Subramanian et al. (2024), Tabel 2 |

---

## 5. Modifikasi & Justifikasi

### 5.1 Penambahan State Glukagon A(t)

**Alasan:** Pada T2D, glukagon memainkan peran patofisiologis yang signifikan.
Sel-α pankreas tidak menekan sekresi glukagon secara normal ketika glukosa naik,
sehingga kadar glukagon tetap tinggi (hyperglucagonemia). Glukagon kemudian menstimulasi
produksi glukosa hepatik (hepatic glucose production/HGP), memperparah hiperglikemia.

Subramanian et al. (2024) membuktikan bahwa parameter `k1` (konstanta supresi glukagon)
secara signifikan lebih rendah pada T2D (0.15 vs 0.25, p = 0.017), dan parameter
`γ2` (amplitudo sekresi glukagon) secara signifikan lebih tinggi (5.0 vs 2.6, p = 0.009).
Tanpa memasukkan `A(t)`, abnormalitas ini tidak dapat direpresentasikan.

**Konsekuensi implementasi:** Mode Normal memiliki 3 state (G, X, I).
Mode T2D memiliki 4 state (G, X, I, A). Matriks Jacobian dan analisis
equilibrium harus diperbarui.

### 5.2 Modifikasi a1: Resistansi Insulin

**Masalah:** Di Mode Normal, sensitivitas insulin terwakili implisit melalui parameter
`p3` (laju aktivasi) dan keseimbangan sistem secara keseluruhan. Di T2D, sensitivitas
insulin perifer turun drastis — sel otot dan adiposa tidak merespons insulin dengan baik.

**Dalam model Subramanian et al. (2024):** Parameter `a1` adalah konstanta laju untuk
insulin-dependent glucose disposal. Pada kontrol sehat, `a1` rata-rata 5.0×10⁻⁴; pada
T2D turun menjadi 1.1×10⁻⁴ (penurunan 4.5 kali lipat, p < 0.0021).

**Solusi:** `a1` di Mode T2D ditetapkan 1.1×10⁻⁴. Ini adalah **parameter kunci** yang
paling membedakan T2D dari kondisi normal dalam model ini.

**Catatan satuan:** Parameter `a1` di Subramanian et al. (2024) dalam satuan
`(10 pM·min)⁻¹`, berbeda dengan konvensi Bergman (mg/dL, μU/mL). Konversi satuan
wajib diperiksa sebelum implementasi (lihat Bagian 9).

### 5.3 Modifikasi k1: Disfungsi α-Cell

**Masalah:** Pada kondisi normal, glukosa yang naik menekan sekresi glukagon secara
eksponensial melalui mekanisme parakrin di islet pankreas. Pada T2D, mekanisme ini
terganggu — glukagon tidak cukup tertekan, sehingga produksi glukosa hepatik berlanjut
bahkan saat glukosa darah sudah tinggi.

**Solusi:** `k1` diturunkan dari 0.25 menjadi 0.15 mM⁻¹. Nilai ini berarti kurva
supresi glukagon lebih datar — glukagon tidak turun secepat pada kondisi normal
meski glukosa naik.

### 5.4 Modifikasi γ3: Penurunan Efek Incretin

**Masalah:** Incretin (GLP-1 dan GIP) adalah hormon usus yang memperkuat sekresi insulin
setelah makan secara oral. Pada T2D, efek incretin secara signifikan melemah — rata-rata
turun dari ~57% (sehat) menjadi ~30% (T2D) dari total respons insulin post-prandial
(Subramanian et al., 2024, Tabel 4 dan 5, p < 0.032).

**Dalam model ini:** `γ3` merepresentasikan laju sekresi insulin yang dipotensiai oleh
incretin. Nilainya turun dari 0.0068 menjadi 0.0017 (mg/dL·min)⁻¹ — penurunan ~4×
yang konsisten dengan data klinis.

**Implikasi simulasi:** Meski respons glukosa terlihat mirip dengan Mode Normal saat puasa,
perbedaan terbesar akan terlihat **setelah makan** — insulin pada T2D naik lebih lambat
dan lebih rendah, menyebabkan hiperglikemia post-prandial lebih lama.

### 5.5 Peningkatan Gb: Fasting Hyperglycemia

**Alasan:** Berbeda dengan T1D (pankreas tidak berfungsi → glukosa tak terbendung),
pada T2D pankreas masih berfungsi sebagian. Equilibrium baru sistem T2D terjadi pada
glukosa basal yang lebih tinggi — umumnya 110–160 mg/dL saat puasa.

Nilai `Gb = 126 mg/dL` dipilih sebagai batas diagnostik T2D (ADA: FPG ≥ 126 mg/dL),
yang dikonfirmasi oleh kondisi awal pada Banzi et al. (2021) (117 ± 7 mg/dL untuk
grup diabetik vs 85 ± 8 mg/dL untuk sehat).

---

## 6. Analisis Equilibrium

Syarat equilibrium untuk sistem 4-state (G, X, I, A) di titik basal
(G=Gb, X=0, I=Ib, A=Ab):

```
Dari persamaan (1) di equilibrium:
  -(SG + a1·Ib)·Gb + a2·Ab + D(t→∞) = 0
  → untuk D=0: (SG + a1·Ib)·Gb = a2·Ab
  → syarat: SG·Gb + a1·Ib·Gb = a2·Ab   ... (*)

Dari persamaan (2):
  -p2·0 + p3·(Ib - Ib) = 0   ✓ otomatis terpenuhi

Dari persamaan (3) di equilibrium:
  -n1·(Ib - Ib) + γ1·ψ(Gb) + γ3·Φ(Gb) = 0
  → γ1·ψ(Gb) + γ3·Gb = n1·0 = 0
  → TIDAK nol kecuali γ1 dan γ3 = 0
  → ini berarti ada net secretion → perlu suku klirens dari I

  Persamaan (3) yang benar di equilibrium:
  0 = -n1·(Ib - Ib) + γ1·ψ(Gb) + γ3·Gb - net_clearance
  → syarat equilibrium: γ1·ψ(Gb) + γ3·Gb = n1·(I_eq - Ib)

Dari persamaan (4) di equilibrium:
  -n2·(Ab - Ab) + γ2·φ(Gb) = 0
  → γ2·φ(Gb) = 0
  → Tidak terpenuhi kecuali φ(Gb) = 0 atau γ2 = 0
```

**Implikasi:** Syarat equilibrium sistem 4-state lebih kompleks dari Mode Normal.
Verifikasi analitik syarat `(*)` dan keseimbangan persamaan insulin dan glukagon
**wajib dilakukan sebelum implementasi** (lihat Bagian 9).

---

## 7. Perilaku Fisiologis yang Diharapkan

Berdasarkan literatur dan hasil simulasi di keempat paper, Mode T2D diharapkan
mereproduksi perilaku berikut yang **berbeda secara kualitatif** dari Mode Normal:

| Perilaku | Mode Normal | Mode T2D yang Diharapkan |
|---|---|---|
| Glukosa basal (puasa) | 92 mg/dL | **126+ mg/dL** |
| Glukosa puncak pasca-makan | ~122 mg/dL | **>160–200 mg/dL** |
| Waktu kembali ke basal | ~180 menit | **>300 menit** (atau tidak kembali) |
| Insulin puncak | ~27 μU/mL | **Lebih rendah atau tertunda** |
| Delay insulin vs glukosa | ~21.5 menit | **Lebih panjang** (fase pertama hilang) |
| Glukagon saat makan | Tertekan normal | **Tidak cukup tertekan** (tetap tinggi) |
| Glukagon basal | ~10–12 pM | **~20 pM** (lebih tinggi) |
| Undershoot setelah makan | Ada (realistis) | Mungkin lebih kecil atau tidak ada |

---

## 8. Perbandingan dengan Keempat Paper Sumber

### 8.1 Subramanian et al. (2024)

Paper ini adalah **sumber utama** parameter dan persamaan T2D dalam model ini.
Sistem persamaan (1)–(3) diadaptasi langsung dari persamaan OGTT paper ini, dengan
penyederhanaan pada fungsi incretin (GLP-1 dan GIP dilebur menjadi suku `γ3`).
Data klinis: 8 pasien T2D, 75 g OGTT, data glukosa–insulin–glukagon 240 menit.
RMSE fitting: R² > 0.98 untuk semua subjek.

**Batas adaptasi:** Paper asli menggunakan delay diferensial (DDE) untuk glukagon
dan incretin. Dalam implementasi ini, delay digantikan dengan parameter tunggal
`τ1` dan diabaikan untuk incretin demi kesederhanaan.

### 8.2 Yang et al. / PLOS Comp. Biol. (2023)

Paper ini memberikan **justifikasi konseptual** bahwa resistansi insulin sebaiknya
dimodelkan sebagai fungsi menurun dari level insulin (`C(I)` dinamis), bukan konstanta.
Dalam model ini, disederhanakan menjadi parameter tetap `a1` yang bernilai rendah untuk
mempertahankan kompatibilitas dengan Mode Normal.

Paper ini juga menunjukkan bahwa model progresi T2D jangka panjang memerlukan state
massa sel-β (`β`) — elemen ini **tidak disertakan** dalam versi awal ini karena scope
simulator adalah simulasi per episode (menit–jam), bukan progresivitas bertahun-tahun.

### 8.3 López-Palau & Olais-Govea (2020)

Paper ini mengidentifikasi secara eksplisit **tiga metabolic rate yang sensitif** dalam
T2D: `M_PGU^I` (sensitivitas insulin perifer), `M_HGP^G/I` (produksi glukosa hepatik),
dan `r_PIR` (pelepasan insulin pankreas). Ketiga abnormalitas ini terpetakan ke
parameter `a1`, `a2·A(t)`, dan `γ1·ψ(G) + γ3·Φ(G)` dalam model ini secara berturutan.

Model 28-ODE paper ini **tidak diadopsi langsung** karena terlalu kompleks untuk
simulator edukasi, namun memberikan validasi fisiologis bahwa tiga jalur tersebut
memang yang paling signifikan di T2D.

### 8.4 Banzi et al. (2021)

Paper ini memberikan **validasi dengan data pasien nyata** (103 pasien, Rwanda) dan
mengkonfirmasi nilai kondisi awal T2D: glukosa jaringan 117 ± 7 mg/dL, insulin
4.5 ± 0.4 μU/mL (vs 85 ± 8 mg/dL dan 5.1 ± 0.3 μU/mL pada sehat).

Konsep **labile insulin** (`Ib`) dan **stored insulin** (`Is`) dari paper ini belum
disertakan dalam model ini. Penambahan dua kompartemen pankreas ini merupakan
kandidat pengembangan iterasi berikutnya.

---

## 9. Konversi Satuan (WAJIB Diperiksa)

Integrasi parameter dari Subramanian et al. (2024) dengan model Bergman memerlukan
konversi satuan yang hati-hati karena kedua sistem menggunakan konvensi berbeda:

| Kuantitas | Satuan Bergman/Mode Normal | Satuan Subramanian (2024) | Catatan |
|---|---|---|---|
| Glukosa G | mg/dL | mM (milimolar) | 1 mM = 18 mg/dL |
| Insulin I | μU/mL | pM (picomolar) | 1 μU/mL ≈ 6.0 pM |
| Glukagon A | — | pM | State baru |
| Parameter a1 | — | (10 pM·min)⁻¹ | Perlu konversi ke satuan mg/dL dan μU/mL |
| Konstanta K (Hill) | — | 17.0 mM = 306 mg/dL | Konversi diperlukan |
| Konstanta k1 | — | mM⁻¹ = (18 mg/dL)⁻¹ | Konversi diperlukan |

**Langkah konversi yang wajib dilakukan sebelum implementasi:**

```
1. Konversi K_Hill: 17.0 mM × 18 = 306 mg/dL (gunakan dalam ψ(G))
2. Konversi a1:
   - a1_paper = 1.1 × 10⁻⁴ (10 pM·min)⁻¹
   - G dalam mg/dL, I dalam μU/mL
   - Perlu faktor konversi: 1 μU/mL ≈ 6.0 pM → 10 pM ≈ 1.667 μU/mL
   - a1_impl = 1.1 × 10⁻⁴ / 1.667 ≈ 6.6 × 10⁻⁵ (μU/mL·min)⁻¹
3. Konversi k1:
   - k1_paper = 0.15 mM⁻¹
   - k1_impl = 0.15 / 18 = 0.00833 (mg/dL)⁻¹
4. Periksa dimensi setiap suku pada persamaan (1)–(4)
   agar semua suku berdimensi mg/(dL·min)
```

---

## 10. Parameter Final untuk Implementasi (Tentatif)

> **⚠️ Status:** nilai-nilai berikut adalah **estimasi awal** berdasarkan konversi
> dari keempat paper. Semua nilai harus diverifikasi secara numerik sebelum digunakan
> dalam simulator UI.

```python
# === Parameter dipertahankan dari Mode Normal ===
SG  = 0.03082      # 1/min (glucose effectiveness, identik p1)
p2  = 0.02093      # 1/min
p3  = 1.062e-5     # L/(min²·mU)
Gb  = 126.0        # mg/dL (fasting T2D, bukan 92 seperti normal)
Ib  = 7.3          # μU/mL (basal insulin; bisa lebih tinggi di early T2D)
G0  = 126.0        # mg/dL
X0  = 0.0
I0  = Ib

# === Parameter T2D baru / dimodifikasi (PERLU KALIBRASI SATUAN) ===
a1  = 6.6e-5       # (μU/mL·min)⁻¹ — RESISTANSI INSULIN (vs 3.0e-4 sehat)
a2  = 0.26         # mg/(dL·pM·min) — efek glukagon hepatik
n1  = 0.14         # 1/min — klirens insulin
n2  = 0.08         # 1/min — klirens glukagon
Ab  = 20.0         # pM — glukagon basal T2D
A0  = Ab

# === Parameter glukagon ===
gamma_2 = 5.0      # pM/min — sekresi glukagon (lebih tinggi dari normal 2.6)
k1_gluc = 0.00833  # (mg/dL)⁻¹ — supresi glukagon (lebih rendah dari normal)
k2_gluc = 0.56     # mM⁻¹ — recovery parameter glukagon

# === Parameter sekresi insulin ===
gamma_1 = 4.0      # 10 pM/min — sekresi insulin bergantung glukosa
gamma_3 = 0.0017   # (mg/dL·min)⁻¹ — INCRETIN EFFECT (vs 0.0068 sehat)

# === Fungsi Hill ===
K_hill  = 306.0    # mg/dL (= 17 mM × 18)
h_hill  = 1.27     # koefisien Hill

# === Input makan (identik Mode Normal) ===
mealTau = 40.0     # menit
# D(t) = A * dtSince * exp(-dtSince/mealTau) / mealTau**2
```

---

## 11. Skenario Demonstrasi

| Skenario | Input | Perilaku yang Diharapkan |
|---|---|---|
| Basal puasa | Tidak ada makan | G tetap ~126 mg/dL, tidak drift |
| Makan normal | Amplitudo makan sedang | Glukosa naik >160 mg/dL, kembali lambat (>3 jam) |
| Makan besar | Amplitudo makan besar | Puncak glukosa bisa >200 mg/dL |
| Bandingkan dengan Normal | Makan identik | T2D puncak lebih tinggi, durasi elevasi lebih panjang |

---

## 12. Verifikasi Wajib Sebelum Implementasi UI

### Lapis 1 — Analitik (Syarat Equilibrium)
Buktikan secara aljabar bahwa titik (G=Gb, X=0, I=Ib, A=Ab) adalah equilibrium
yang valid untuk sistem 4-ODE dengan nilai parameter yang dipilih. Jika tidak valid,
sesuaikan nilai parameter atau tambahkan suku offset.

### Lapis 2 — Numerik: Uji Diam
Jalankan simulasi 1000 menit dari kondisi awal tanpa gangguan makan. Seluruh state
tidak boleh drift. Kriteria lulus: deviasi maksimum dari baseline < 0.01 mg/dL untuk G.

### Lapis 3 — Numerik: Uji Gangguan Makan
Satu episode makan, simulasi 1000 menit. Verifikasi:
- Glukosa puncak > Mode Normal pada amplitudo makan sama
- Durasi elevasi glukosa > Mode Normal
- Glukagon tidak turun normal (tetap lebih tinggi dari Ab)
- Sistem kembali (atau mendekati) baseline setelah waktu yang lama

### Lapis 4 — Perbandingan Kualitatif vs Literatur
Bandingkan profil glukosa–insulin–glukagon simulasi dengan Figure 5 Subramanian et al.
(2024) — data pasien T2D 1 dan 8. Ini adalah validasi kualitatif, bukan fitting.
RMSE yang diharapkan: < 30 mg/dL (tidak fitting satu pasien, hanya verifikasi bentuk).

### Lapis 5 — Uji Konvergensi Langkah Integrasi
Jalankan dengan dt = 0.5, 0.25, dan 0.1 menit. Puncak glukosa tidak boleh berbeda
lebih dari 0.01 mg/dL antar resolusi.

### Ringkasan Status Verifikasi (akan diisi)

| # | Uji | Target | Status |
|---|---|---|---|
| 1 | Analitik (syarat equilibrium) | Terbukti aljabar | 🔲 Belum |
| 2 | Uji diam 1000 menit | Deviasi < 0.01 mg/dL | 🔲 Belum |
| 3 | Uji gangguan makan | G_peak > 160 mg/dL, durasi > 3 jam | 🔲 Belum |
| 4 | Validasi kualitatif vs Subramanian | RMSE < 30 mg/dL | 🔲 Belum |
| 5 | Konvergensi langkah integrasi | Δ < 0.01 mg/dL | 🔲 Belum |
| 6 | Konversi satuan eksplisit | Semua suku berdimensi benar | 🔲 Belum |

---

## 13. Batas Model dan Catatan Penting

1. **Bukan representasi ketoasidosis.** Model ini tidak mencakup metabolisme asam lemak
   bebas, badan keton, atau komplikasi akut. Tidak untuk estimasi risiko klinis.

2. **Progresi T2D tidak dimodelkan.** Model ini adalah *snapshot* kondisi T2D yang sudah
   terbentuk. Progresi bertahun-tahun (seperti pada Yang et al., 2023) memerlukan
   penambahan state massa sel-β (`β`) yang tidak ada di sini.

3. **Incretin disederhanakan.** GLP-1 dan GIP dilebur menjadi satu suku `γ3·Φ(G)`.
   Efek GIP yang merangsang glukagon di kadar glukosa rendah tidak dimodelkan.

4. **Glukagon menggunakan model hiperbolik sederhana.** Perilaku histeresis glukagon
   yang diamati Subramanian et al. (2024) — kurva supresi berbeda dari kurva pemulihan —
   digantikan dengan fungsi eksponensial tunggal `φ(G) = exp(-k1·G)`. Ini adalah
   penyederhanaan yang dapat diperbaiki di iterasi berikutnya.

5. **Parameter belum dikalibrasi sebagai sistem terintegrasi.** Nilai `a1`, `a2`,
   `γ1`, `γ2`, `γ3` masing-masing diambil dari konteks paper yang berbeda dan belum
   dioptimasi bersama-sama. Kalibrasi simultan terhadap satu dataset T2D wajib
   dilakukan sebelum model dinyatakan valid.

---

*Dokumen ini merupakan spesifikasi riset Mode T2D versi awal. Verifikasi penuh
dan kalibrasi parameter adalah prasyarat sebelum integrasi ke simulator UI.
Artificial Pancreas / kontroler PID untuk T2D adalah tahap terpisah sesudah
model T2D open-loop ini selesai diverifikasi.*
