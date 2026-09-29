# Memahami Simulator Glukosa–Insulin dan Hubungannya dengan Dunia Nyata

Dokumen ini menjelaskan makna biologis dan klinis dari komponen yang terlihat pada
simulator, mulai dari tubuh sehat, diabetes tipe 1 dengan terapi manual, hingga
artificial pancreas dengan pengendali PID-IFB. Tujuannya agar pengguna memahami
apa yang sedang divisualisasikan tanpa harus lebih dahulu membaca persamaan
matematis lengkap.

> **Batas penggunaan:** simulator ini merupakan media pendidikan. Subjek yang
> disimulasikan adalah subjek virtual berbasis parameter populasi. Nilai glukosa,
> insulin, makanan, bolus, dan keluaran pengendali bukan prediksi untuk individu,
> bukan alat diagnosis, dan bukan rekomendasi dosis insulin.

---

## 1. Gambaran Besar

Tubuh membutuhkan glukosa sebagai sumber energi. Setelah makanan dicerna,
karbohidrat diubah menjadi glukosa dan masuk ke aliran darah. Insulin membantu
sel, terutama jaringan otot dan lemak, mengambil dan menggunakan glukosa. Hati
juga menyimpan atau melepaskan glukosa sesuai kondisi tubuh.

Hubungan dasarnya dapat dibaca seperti ini:

```text
makanan → glukosa darah naik → insulin bekerja → jaringan mengambil glukosa
                                      ↓
                              glukosa turun kembali
```

Perbedaan utama ketiga mode simulator adalah pihak yang menyediakan dan mengatur
insulin:

| Mode | Sumber insulin dalam simulator | Cara pengaturan |
|---|---|---|
| Normal | Pankreas virtual | Otomatis sebagai respons fisiologis terhadap glukosa |
| Diabetes Tipe 1 | Insulin basal dan bolus dari luar tubuh | Basal tetap, bolus diberikan pengguna |
| Diabetes Tipe 1 + AP | Pompa insulin virtual | Laju pompa diubah otomatis oleh PID-IFB berdasarkan sensor glukosa |

---

## 2. Istilah Dasar yang Muncul di Web

### 2.1 Glukosa

Kartu **Glukosa** menampilkan konsentrasi glukosa darah dalam `mg/dL`. Angka ini
adalah keluaran model pada waktu simulasi saat ini. Setelah tombol makan ditekan,
glukosa tidak muncul seluruhnya sekaligus; simulator memodelkan penyerapan makanan
secara bertahap.

### 2.2 Insulin

Kartu **Insulin** menampilkan konsentrasi insulin plasma dalam `µU/mL`. Ini berbeda
dari dosis insulin dalam `U` yang dimasukkan pengguna:

- `U` adalah jumlah insulin yang diberikan;
- `µU/mL` adalah konsentrasi insulin yang telah mencapai plasma;
- `U/jam` adalah laju pemberian insulin oleh pompa.

Ketiganya tidak dapat disamakan langsung karena insulin yang diberikan harus
diserap dahulu dan tersebar dalam tubuh.

### 2.3 Efek insulin di jaringan

Insulin plasma tidak langsung menurunkan glukosa pada detik yang sama. Model
memiliki variabel aksi insulin yang menggambarkan keterlambatan antara insulin
masuk, beredar, lalu meningkatkan pemakaian glukosa dan menekan produksi glukosa
oleh hati. Karena itu, kurva insulin dan glukosa tidak bergerak bersamaan.

### 2.4 Makan

Pada mode T1D dan AP, preset makanan adalah:

| Pilihan | Karbohidrat dalam simulator |
|---|---:|
| Kecil | 25 g |
| Sedang | 50 g |
| Besar | 75 g |

Karbohidrat tersebut merupakan gangguan yang menaikkan laju kemunculan glukosa.
Model juga memakai fraksi absorpsi dan waktu penyerapan. Jadi, `50 g` karbohidrat
tidak berarti 50 gram glukosa langsung dimasukkan ke darah.

Pada **mode Normal**, preset kecil/sedang/besar masih memakai amplitudo gangguan
model, bukan konversi gram karbohidrat. Karena itu, perbandingan kuantitatif ukuran
makanan antara mode Normal dan mode T1D/AP tidak boleh dianggap setara.

### 2.5 Status glikemik

Status pada web mengelompokkan nilai glukosa saat ini agar perubahan mudah
diamati. Status tersebut adalah penanda visual sederhana, bukan diagnosis medis.

---

## 3. Mode Normal: Pankreas Sehat

### 3.1 Apa yang terjadi di dunia nyata?

Pada orang tanpa diabetes, pankreas memiliki sel beta yang melepaskan insulin.
Ketika glukosa naik setelah makan, sekresi insulin bertambah. Insulin membantu
jaringan mengambil glukosa dan membantu mengendalikan pelepasan glukosa oleh
hati. Saat glukosa kembali turun, sekresi insulin juga berkurang.

Artinya, tubuh memiliki umpan balik alami:

```text
glukosa diukur oleh sistem biologis
        ↓
pankreas menyesuaikan insulin
        ↓
insulin memengaruhi glukosa
        └──────── umpan balik ────────┘
```

### 3.2 Apa yang dilakukan simulator?

Mode Normal memakai model glukosa–insulin tiga-state berbasis Bergman–Pacini yang
telah digabungkan menjadi sistem loop tertutup untuk keperluan pendidikan:

- `G(t)` menggambarkan glukosa plasma;
- `I(t)` menggambarkan insulin plasma;
- `X(t)` menggambarkan efek insulin aktif di jaringan;
- `D(t)` menggambarkan gangguan makanan;
- pankreas virtual menaikkan insulin ketika glukosa melewati ambang respons.

Simulator mulai pada keadaan basal sekitar `92 mg/dL` untuk glukosa dan
`7,3 µU/mL` untuk insulin. Angka ini adalah parameter model, bukan nilai target
untuk semua orang.

### 3.3 Cara mencoba

1. Pilih mode **Normal**.
2. Tekan **Makan** dan pilih ukuran makanan.
3. Amati glukosa naik, kemudian insulin pankreas virtual ikut meningkat.
4. Amati efek insulin `X(t)` meningkat setelah insulin plasma.
5. Biarkan simulasi berjalan untuk melihat sistem kembali mendekati baseline.

Pelajaran utamanya adalah tubuh sehat mengatur insulin secara internal tanpa
pengguna memasukkan dosis.

---

## 4. Diabetes Tipe 1: Terapi Basal–Bolus Manual

### 4.1 Apa yang berubah pada diabetes tipe 1?

Pada diabetes tipe 1, kerusakan autoimun pada sel beta menyebabkan produksi
insulin endogen sangat sedikit atau tidak mencukupi. Karena insulin tidak lagi
tersedia secara memadai, insulin perlu diberikan dari luar tubuh.

Terapi sehari-hari biasanya memiliki dua fungsi:

- **basal:** memenuhi kebutuhan insulin dasar di antara waktu makan dan saat tidur;
- **bolus:** tambahan insulin untuk makanan atau koreksi sesuai rencana terapi.

### 4.2 Basal di dunia nyata

Kebutuhan basal dapat dipenuhi dengan dua pendekatan umum:

1. **Suntikan insulin kerja panjang.** Insulin disuntikkan pada jadwal tertentu dan
   diserap perlahan selama berjam-jam.
2. **Pompa insulin.** Pompa memberikan insulin kerja cepat dalam jumlah kecil
   secara kontinu, biasanya dinyatakan dalam `U/jam`.

Checkbox **Infus basal aktif** di simulator paling dekat dengan pendekatan pompa:
selama aktif, model menerima laju insulin dasar yang konstan. Simulator belum
memiliki model farmakokinetik khusus untuk sediaan insulin kerja panjang. Suntikan
kerja panjang dapat memenuhi fungsi basal di dunia nyata, tetapi tidak boleh
dianggap identik secara matematis dengan infus konstan ini.

Basal pada mode ini bukan PID dan tidak menyesuaikan diri terhadap glukosa:

```text
basal aktif    → laju basal tetap diberikan
basal nonaktif → laju basal menjadi nol
```

Mematikan basal tersedia sebagai eksperimen edukasi untuk memperlihatkan akibat
hilangnya insulin dasar. Ini bukan petunjuk untuk menghentikan terapi nyata.

### 4.3 Bolus di dunia nyata

Bolus biasanya memakai insulin kerja cepat. Dalam praktik klinis, dosis dapat
dipengaruhi oleh jumlah karbohidrat, glukosa saat ini, sensitivitas insulin,
insulin yang masih aktif, aktivitas fisik, kondisi sakit, dan instruksi tenaga
kesehatan. Simulator tidak melakukan perhitungan dosis personal tersebut.

Input **2 U** di web berarti pengguna memberikan total 2 unit insulin kepada
subjek virtual. Karena model Hovorka dinormalisasi terhadap berat badan, simulator
mengubahnya secara internal:

```text
2 U = 2.000 mU
berat subjek virtual = 70 kg
2.000 mU / 70 kg = 28,57 mU/kg
```

Massa `28,57 mU/kg` ditambahkan satu kali ke depot insulin subkutan pertama
`S1`. Insulin kemudian berpindah melalui `S1 → S2 → plasma` dan baru menghasilkan
aksi insulin. Jadi, 2 U tidak langsung menjadi konsentrasi plasma dan tidak
langsung menurunkan glukosa pada saat tombol ditekan.

Angka 70 kg adalah asumsi subjek virtual dari konfigurasi simulator. Angka 2 U
adalah preset eksperimen untuk skenario 50 g pada model ini, bukan rekomendasi
dosis untuk orang dengan berat 70 kg.

### 4.4 Apa yang dilakukan model Hovorka?

Mode T1D memakai delapan state utama:

| State | Makna sederhana |
|---|---|
| `Q1` | Massa glukosa pada kompartemen yang dapat diakses |
| `Q2` | Massa glukosa pada kompartemen perifer |
| `S1`, `S2` | Dua tahap penyerapan insulin subkutan |
| `I` | Insulin plasma |
| `x1` | Efek insulin pada distribusi/pemakaian glukosa |
| `x2` | Efek insulin pada pemakaian glukosa perifer |
| `x3` | Efek insulin pada produksi glukosa endogen |

Makanan menaikkan kemunculan glukosa, sedangkan basal dan bolus masuk melalui
depot insulin. Model menghitung perubahan massa glukosa, insulin plasma, dan aksi
insulin dari waktu ke waktu.

### 4.5 Cara mencoba

1. Pilih **Diabetes**.
2. Biarkan **Infus basal aktif**.
3. Pilih makan sedang, yaitu 50 g karbohidrat.
4. Amati kenaikan glukosa tanpa bolus selama beberapa saat.
5. Masukkan `2` pada bolus lalu tekan **Beri Bolus**.
6. Amati insulin diserap bertahap dan glukosa kemudian dipengaruhi oleh aksi insulin.
7. Tekan **Reset**, ulangi tanpa bolus, lalu bandingkan kurvanya.
8. Untuk eksperimen terpisah, matikan basal dan amati perubahan jangka panjang.

Pelajaran utamanya adalah pada T1D insulin berasal dari luar tubuh, dan waktu serta
jumlah insulin memengaruhi respons model.

---

## 5. Diabetes Tipe 1 + PID-IFB: Artificial Pancreas

### 5.1 Apa itu artificial pancreas?

Artificial pancreas atau sistem closed-loop menghubungkan tiga komponen:

```text
sensor glukosa kontinu (CGM) → algoritma pengendali → pompa insulin
             ↑                                      ↓
             └──────── respons tubuh pengguna ──────┘
```

CGM mengukur glukosa secara berkala. Algoritma memakai hasil sensor untuk
menentukan laju insulin, lalu pompa mengirim insulin. Perubahan glukosa berikutnya
dibaca kembali oleh sensor sehingga terbentuk loop tertutup.

Sistem nyata memiliki lebih banyak pengaman, batas perangkat, kalibrasi, alarm,
dan validasi daripada simulator ini.

### 5.2 Sensor CGM dalam simulator

Sensor tidak dianggap membaca glukosa plasma secara instan. Simulator memakai
sensor orde pertama dengan keterlambatan sehingga nilai sensor bergerak menuju
nilai glukosa darah secara bertahap. Ini menggambarkan secara sederhana adanya
lag antara glukosa darah dan pembacaan sensor.

Simulator saat ini belum mensimulasikan seluruh gangguan perangkat nyata seperti
noise, bias, kehilangan data, kesalahan kalibrasi, atau kegagalan set infus.

### 5.3 Apa arti PID?

Pengendali membandingkan glukosa sensor `SG` dengan target `120 mg/dL`:

```text
error = SG − target
```

PID memiliki tiga komponen:

- **P — proportional:** bereaksi terhadap selisih glukosa saat ini;
- **I — integral:** mengakumulasi selisih yang bertahan dari waktu ke waktu;
- **D — derivative:** bereaksi terhadap arah dan kecepatan perubahan glukosa.

Dalam simulator, pengendali diperbarui setiap satu menit. Keluaran akhirnya adalah
laju pompa dalam `U/jam`, bukan bolus manual.

### 5.4 Apa arti IFB?

IFB adalah **insulin feedback**. Pengendali memperkirakan insulin yang masih aktif
dari riwayat insulin yang telah diberikan. Estimasi ini digunakan sebagai umpan
balik agar pengendali tidak terus menambah insulin seolah-olah dosis sebelumnya
sudah tidak ada.

Secara konsep:

```text
PID meminta insulin berdasarkan error glukosa
IFB memperhitungkan insulin yang diperkirakan masih bekerja
hasil gabungan menentukan laju pompa
```

IFB pada simulator adalah estimator matematis tiga-kompartemen, bukan pengukuran
langsung kadar insulin pasien.

### 5.5 Pengaman yang dimodelkan

Implementasi memiliki beberapa pembatas sederhana:

- laju pompa tidak boleh negatif;
- laju pompa dibatasi maksimum;
- integral dibatasi dan dapat dibekukan ketika keluaran jenuh;
- insulin dihentikan sementara jika pembacaan sensor berada pada atau di bawah
  `70 mg/dL`.

Fitur tersebut membantu kestabilan simulasi, tetapi tidak cukup untuk menjadikan
algoritma ini perangkat medis.

### 5.6 Cara mencoba

1. Pilih **Diabetes + AP**.
2. Perhatikan target, basal, laju pompa, dan status pengendali.
3. Pilih makanan. Pengendali simulator tidak diberi tahu bahwa makanan terjadi.
4. Amati glukosa naik dan laju pompa berubah setelah sensor menangkap perubahan.
5. Buka **Perhitungan Matematis** untuk melihat error, komponen P, I, D, IFB, dan
   laju pompa.
6. Bandingkan hasilnya dengan mode Diabetes manual tanpa bolus dan dengan bolus.

Mode AP ini bersifat **fully closed-loop tanpa informasi makan** dalam konteks
simulator. Karena tidak ada pengumuman makanan atau bolus pramakan, lonjakan awal
glukosa dapat tetap besar sebelum pengendali bereaksi.

---

## 6. Perbandingan Tiga Mode

| Pertanyaan | Normal | T1D manual | T1D + PID-IFB |
|---|---|---|---|
| Apakah pankreas virtual menghasilkan insulin? | Ya | Tidak | Tidak |
| Dari mana insulin berasal? | Sekresi pankreas model | Basal tetap + bolus pengguna | Pompa yang diatur controller |
| Siapa yang bereaksi terhadap makanan? | Pankreas model | Pengguna melalui bolus | Controller bereaksi setelah sensor berubah |
| Apakah ada sensor CGM? | Tidak | Tidak | Ya, model sensor dengan lag |
| Apakah laju insulin berubah otomatis? | Sekresi berubah otomatis | Basal tidak berubah otomatis | Ya, setiap satu menit |
| Model utama | Bergman–Pacini termodifikasi | Hovorka delapan-state | Hovorka + sensor + PID-IFB |

---

## 7. Eksperimen Edukasi yang Disarankan

Lakukan setiap eksperimen setelah menekan **Reset** agar kondisi awal sama.

### Eksperimen A — Respons pankreas sehat

Jalankan satu makanan di mode Normal. Perhatikan urutan glukosa naik, insulin
naik, efek insulin meningkat, lalu glukosa kembali mendekati baseline.

### Eksperimen B — T1D tanpa bolus

Pada mode Diabetes dengan basal aktif, berikan makan 50 g tanpa bolus. Bandingkan
puncak dan lamanya glukosa tinggi dengan mode Normal. Perbandingan bentuk kurva
berguna secara konseptual, sedangkan ukuran makan antarmode tidak setara secara
kuantitatif.

### Eksperimen C — T1D dengan bolus

Ulangi makanan 50 g dan berikan bolus 2 U. Bandingkan dengan eksperimen B. Fokus
pada keterlambatan antara pemberian bolus, kenaikan insulin plasma, dan perubahan
glukosa.

### Eksperimen D — Basal dihentikan

Reset mode Diabetes, matikan basal, dan jangan beri makan. Amati insulin serta aksi
insulin berkurang dan glukosa naik. Eksperimen ini menjelaskan fungsi insulin dasar,
bukan tindakan yang boleh dicoba pada terapi nyata.

### Eksperimen E — Closed-loop tanpa informasi makan

Pada mode AP, berikan makan 50 g. Amati sensor, respons laju pompa, dan glukosa.
Bandingkan dengan T1D manual. Ini menunjukkan manfaat dan keterbatasan pengendali
yang baru bereaksi setelah perubahan glukosa terdeteksi.

---

## 8. Hal yang Sengaja Disederhanakan

Simulator belum merepresentasikan seluruh fisiologi atau perawatan diabetes.
Beberapa penyederhanaan penting adalah:

- satu subjek virtual dengan berat tetap 70 kg;
- parameter populasi, bukan parameter hasil identifikasi individu;
- tidak ada glukagon dan respons kontraregulasi hipoglikemia yang lengkap;
- tidak ada aktivitas fisik, stres, sakit, tidur, hormon menstruasi, atau variasi harian;
- tidak ada pemodelan khusus insulin kerja panjang;
- tidak ada komposisi makanan selain jumlah karbohidrat pada mode T1D/AP;
- tidak ada ketidakpastian penghitungan karbohidrat;
- model CGM belum mencakup noise, bias, dropout, dan kegagalan sensor;
- model pompa belum mencakup sumbatan, kebocoran, atau kegagalan set infus;
- PID-IFB merupakan rancangan edukasi yang belum divalidasi sebagai controller
  untuk penggunaan klinis.

Karena itu, kurva harus dibaca sebagai demonstrasi hubungan sebab-akibat dalam
model, bukan prediksi pasti tubuh manusia.

---

## 9. Sumber Ilmiah Utama

1. Bergman, R.N., Ider, Y.Z., Bowden, C.R., & Cobelli, C. (1979).
   *Quantitative estimation of insulin sensitivity*. American Journal of
   Physiology, 236(6), E667–E677.
2. Pacini, G., & Bergman, R.N. (1986). *MINMOD: a computer program to calculate
   insulin sensitivity and pancreatic responsivity from the frequently sampled
   intravenous glucose tolerance test*. Computer Methods and Programs in
   Biomedicine, 23(2), 113–122.
3. Hovorka, R., Canonico, V., Chassin, L.J., et al. (2004). *Nonlinear model
   predictive control of glucose concentration in subjects with type 1 diabetes*.
   Physiological Measurement, 25(4), 905–920.
   DOI: `10.1088/0967-3334/25/4/010`.
4. Steil, G.M., Rebrin, K., Darwin, C., Hariri, F., & Saad, M.F. (2006).
   *Feasibility of automating insulin delivery for the treatment of type 1
   diabetes*. Diabetes, 55(12), 3344–3350.
5. Steil, G.M., Palerm, C.C., Kurtz, N., et al. (2011). *The effect of insulin
   feedback on closed loop glucose control*. Journal of Clinical Endocrinology &
   Metabolism, 96(5), 1402–1408.
6. Ruiz, J.L., Sherr, J.L., Cengiz, E., et al. (2012). *Effect of insulin feedback
   on closed-loop glucose control: a crossover study*. Journal of Diabetes Science
   and Technology, 6(5), 1123–1130.
   DOI: `10.1177/193229681200600517`.

Persamaan, parameter, modifikasi proyek, dan hasil verifikasi numerik dijelaskan
lebih rinci dalam `model-normal-glucose-insulin.md` pada folder yang sama.

