# Rencana Implementasi Panel Parameter untuk Semua Mode

## Tujuan

Panel `Parameter Fisiologis` yang saat ini terutama tampil pada Mode Normal perlu tersedia pada semua mode simulator:

- Normal
- Diabetes Tipe 1
- Diabetes Tipe 1 + AP
- Diabetes Tipe 2
- T2D + PI

Pengguna dapat mengubah parameter yang memiliki rentang numerik tervalidasi. Setelah perubahan diterapkan, simulator menghitung ulang kondisi awal dan menyinkronkan seluruh UI. Persamaan model, solver, dan controller tetap menjadi sumber kebenaran; pekerjaan ini menambah lapisan konfigurasi dan UI.

## Batasan implementasi

1. Jangan mengganti persamaan model atau algoritma controller sebagai bagian dari pekerjaan panel.
2. Jangan mencampur state dari mode sebelumnya dengan mode yang baru dipilih.
3. Jangan menerapkan parameter baru ke simulasi yang sedang berjalan secara diam-diam.
4. Parameter yang belum memiliki rentang aman harus ditampilkan sebagai read-only atau ditunda.
5. Nilai parameter dan deskripsi harus menyebut satuan serta arti edukatifnya.

## Struktur data parameter

Buat registry UI terpisah dari konstanta model. Setiap definisi parameter minimal memiliki:

```js
{
  key: 'p3',
  symbol: 'p₃',
  label: 'Sensitivitas insulin',
  unit: 'per menit per µU/mL',
  default: 1.062e-5,
  min: 0.3e-5,
  max: 3.0e-5,
  step: 0.02e-5,
  decimals: 6,
  group: 'physiology',
  editable: true,
  description: 'Mengatur seberapa kuat insulin menurunkan glukosa.'
}
```

Registry harus memiliki definisi terpisah untuk setiap mode. Jangan menggunakan satu daftar parameter untuk semua mode karena state dan satuannya berbeda.

Simpan nilai edit sementara di objek UI/model terpisah dari `simState`. Nilai baru baru diterapkan setelah pengguna menekan tombol `Terapkan & Reset`.

## Parameter per mode

### Mode Normal

Parameter utama yang dapat diedit:

- `p1` — efektivitas glukosa
- `p3` — sensitivitas insulin
- `p6` — respons pankreas

Parameter lanjutan:

- `p2` — peluruhan efek insulin
- `n` — pembersihan insulin
- `p5` — ambang sekresi pankreas
- `Gb` — glukosa basal
- `Ib` — insulin basal
- `mealTau` — waktu penyerapan makanan

### Diabetes Tipe 1

Parameter utama:

- `weightKg` — berat badan model
- `tmaxG` — waktu penyerapan glukosa dari makanan
- `tmaxI` — waktu penyerapan insulin subkutan

Parameter lanjutan:

- `SIT` — sensitivitas transport glukosa
- `SID` — sensitivitas penggunaan glukosa
- `SIE` — sensitivitas produksi glukosa hati
- `EGP0` — produksi glukosa endogen
- `F01` — penggunaan glukosa non-insulin
- `VG` — volume distribusi glukosa
- `VI` — volume distribusi insulin
- `k12` — perpindahan glukosa antarkompartemen
- `ka1`, `ka2`, `ka3` — laju penyerapan insulin

Kontrol terapi yang tetap berada di panel skenario:

- bolus insulin
- toggle basal aktif

### Diabetes Tipe 1 + AP

Gunakan parameter fisiologis T1D yang sama. Tambahkan kelompok `controller`:

Parameter utama:

- target glukosa
- `Kp`
- `Ti`
- `Td`

Parameter lanjutan:

- laju pompa maksimum
- batas low-glucose suspend
- konstanta waktu sensor
- periode sampling controller
- parameter insulin feedback

### Diabetes Tipe 2

Parameter utama:

- skala resistensi insulin
- `p3` dasar
- `p6` respons sel beta
- `Gb` glukosa basal

Parameter lanjutan:

- `p1`
- `p2`
- `n`
- `p5`
- `Ib`
- `mealTau`

Model menggunakan:

```text
p3_efektif = p3_dasar / resistanceScale
```

UI harus menampilkan `p3 efektif` setelah slider resistensi berubah agar dampaknya dapat dipahami.

### T2D + PI

Gunakan parameter fisiologis T2D yang sama. Tambahkan kelompok `controller`:

Parameter utama:

- target glukosa
- `Kp`
- `Ti`

Parameter lanjutan:

- laju insulin maksimum
- batas suspend
- konstanta waktu sensor
- `tmaxI` aktuator subkutan
- berat badan aktuator
- periode sampling PI

Bias controller harus dihitung ulang dari parameter plant, target, berat badan, dan skala resistensi.

## Layout UI

Panel mengikuti pola yang sudah digunakan pada Mode Normal:

```text
Parameter Model
├─ Parameter fisiologis
│  ├─ 3–5 parameter utama
│  └─ Parameter lainnya (details)
├─ Parameter controller (AP/PI saja)
│  ├─ parameter utama
│  └─ Parameter controller lainnya (details)
└─ [Default Mode] [Terapkan & Reset]
```

Setiap baris parameter menampilkan:

- simbol dan nama;
- nilai saat ini;
- nilai default;
- satuan;
- slider atau input angka;
- deskripsi singkat;
- penanda `diubah` jika nilainya berbeda dari default.

Pada layar ponsel, kelompok menjadi satu kolom. Tombol penerapan tetap terlihat setelah daftar parameter.

## Alur penerapan perubahan

1. Pengguna menggeser slider atau mengisi nilai.
2. Nilai divalidasi terhadap `min`, `max`, dan `step`.
3. UI menandai bahwa ada perubahan yang belum diterapkan.
4. Tombol `Terapkan & Reset` menjadi aktif.
5. Sistem menyalin nilai edit ke konfigurasi mode aktif.
6. Sistem menghitung ulang equilibrium dengan fungsi model yang sudah ada.
7. Jika equilibrium valid, panggil reset mode secara penuh.
8. Jika gagal, batalkan penerapan, pertahankan simulasi sebelumnya, dan tampilkan pesan parameter yang bermasalah.

Reset penuh harus menyegarkan sekaligus:

- `simState`;
- state controller;
- basal dan bias;
- grafik utama serta grafik matematika;
- status KPI;
- label target;
- panel anatomi dan particle flow;
- ringkasan grafik;
- nilai input dan slider.

## Penyimpanan per mode

Setiap mode memiliki draft parameter dan nilai aktif sendiri. Ketika pengguna berganti mode:

1. simpan draft mode sebelumnya;
2. pilih registry mode baru;
3. tampilkan nilai aktif mode baru;
4. reset simulator mode baru;
5. perbarui subtitle, panel kontrol, grafik, dan persamaan.

Mengubah parameter T2D tidak boleh mengubah parameter Normal, T1D, AP, atau PI.

## Validasi dan batasan

Validasi wajib mencakup:

- angka finite;
- nilai tidak boleh negatif jika state/model mensyaratkannya;
- batas minimum dan maksimum;
- kelipatan step;
- target controller berada pada rentang glukosa yang masuk akal untuk model;
- `Ti`, waktu penyerapan, dan konstanta waktu tidak boleh nol;
- laju pompa tidak boleh melebihi batas controller.

Parameter yang dapat membuat equilibrium tidak stabil sebaiknya dimulai sebagai read-only. Setelah rentangnya diuji, parameter tersebut baru dapat diaktifkan.

## Sinkronisasi grafik dan anotasi

- Mode AP menampilkan garis target AP pada grafik.
- Mode T2D + PI menampilkan garis target PI.
- Perubahan target langsung memperbarui label garis dan panel controller setelah diterapkan.
- Sumbu grafik dihitung ulang berdasarkan parameter dan riwayat baru.
- Panel matematika menampilkan nilai parameter aktif, bukan nilai default.

## Sinkronisasi anatomi

- Normal dan T2D menampilkan insulin endogen dari pankreas.
- T1D dan AP menampilkan insulin eksogen dari depot subkutan.
- T2D + PI dapat menampilkan insulin dari pankreas dan depot subkutan.
- Deskripsi aksesibilitas SVG harus mengikuti mode aktif.

## Checklist implementasi

- [ ] Buat registry parameter per mode.
- [ ] Pisahkan draft parameter dari `simState`.
- [ ] Render panel parameter berdasarkan mode aktif.
- [ ] Tambahkan kelompok `Parameter lainnya`.
- [ ] Tambahkan kelompok controller untuk AP dan PI.
- [ ] Tambahkan validasi min, max, step, dan satuan.
- [ ] Tambahkan status perubahan yang belum diterapkan.
- [ ] Implementasikan `Terapkan & Reset`.
- [ ] Implementasikan reset parameter mode aktif.
- [ ] Pastikan reset mengosongkan grafik dan riwayat lama.
- [ ] Pastikan bias AP/PI dihitung ulang.
- [ ] Pastikan panel matematika memakai parameter aktif.
- [ ] Pastikan label makan dan anatomi mengikuti mode.
- [ ] Periksa tampilan desktop dan lebar 390 px.
- [ ] Periksa keyboard focus dan reduced motion.

## Checklist verifikasi perilaku

Untuk setiap mode, periksa:

1. buka panel parameter;
2. ubah satu parameter utama;
3. pastikan nilai draft berubah dan tombol penerapan aktif;
4. terapkan perubahan;
5. pastikan simulator kembali ke state awal yang baru;
6. pastikan grafik, KPI, controller, persamaan, dan anatomi konsisten;
7. reset ke default;
8. ganti mode dan pastikan parameter mode lain tidak ikut berubah.

Untuk AP dan PI, tambahan pemeriksaan:

- target baru muncul pada garis grafik;
- bias controller dihitung ulang;
- laju awal tidak memakai nilai controller dari mode sebelumnya;
- suspend dan saturasi tetap memiliki batas;
- perubahan resistensi tidak menghasilkan state negatif.

## Status dan konteks edukasi

Panel ini untuk eksperimen edukasi terhadap model matematika. Parameter, target, bias, dan keluaran controller tidak boleh diberi label sebagai rekomendasi klinis atau dosis pasien nyata.

