# UCSF ChimeraX - MD Trajectory Visualization Guide

Panduan lengkap untuk memvisualisasikan hasil simulasi Molecular Dynamics (MD) menggunakan UCSF ChimeraX dan menghasilkan video Full HD untuk publikasi.

---

## Prasyarat

### File yang Dibutuhkan

1. **Struktur** (salah satu):
   - `md.gro` (GROMACS format)
   - `md.pdb` (PDB format)

2. **Trajectory** (sudah di-downsample):
   - `visual_ready_*.xtc` (hasil proses dari GROMACS)

### Software

- UCSF ChimeraX (Download: <https://www.cgl.ucsf.edu/chimerax/download.html>)

---

## Tahap 1: Persiapan Trajectory di WSL (Wajib!)

Sebelum buka di ChimeraX, trajectory harus "dikecilkan" dulu agar tidak crash.

### Install GROMACS (Sekali Saja)

```bash
sudo apt update
sudo apt install gromacs
```

### Proses Trajectory

```bash
# Masuk ke folder project
cd /mnt/e/Project/herbal-network-pharmacology

# === Untuk 264THM-PPARG ===
cd kaggle-output-sim/264THM_PPARG/md/

# Downsample (skip tiap 10 frame) + Fit protein ke tengah
echo "1 0" | gmx trjconv -s md.tpr -f full_trajectory.xtc \
    -o visual_ready_264THM.xtc -skip 10 -fit rot+trans

# Fix PBC (supaya ligan tidak loncat-loncat)
echo "1 0" | gmx trjconv -s md.tpr -f visual_ready_264THM.xtc \
    -o visual_ready_264THM_fixed.xtc -pbc mol -center

cd /mnt/e/Project/herbal-network-pharmacology

# === Untuk Luteolin-PDE5A ===
cd kaggle-output-luteolin/Luteolin_PDE5A/md/

echo "1 0" | gmx trjconv -s md.tpr -f full_trajectory.xtc \
    -o visual_ready_Luteolin.xtc -skip 10 -fit rot+trans

echo "1 0" | gmx trjconv -s md.tpr -f visual_ready_Luteolin.xtc \
    -o visual_ready_Luteolin_fixed.xtc -pbc mol -center

cd /mnt/e/Project/herbal-network-pharmacology
```

**Output:**

- `visual_ready_264THM_fixed.xtc`
- `visual_ready_Luteolin_fixed.xtc`

---

## Tahap 2: Buka File di ChimeraX

### Langkah

1. Buka **UCSF ChimeraX**
2. **File > Open** → Pilih `md.gro` (struktur)
3. **File > Open** → Pilih `visual_ready_*_fixed.xtc` (trajectory)
4. Slider "Coordinate sets" akan muncul di bawah

---

## Tahap 3: Styling (Tampilan Profesional)

Ketik perintah ini satu per satu di **Command Line** (bawah layar):

```chimerax
# Background putih (standar publikasi)
set bg white

# Sembunyikan air dan ion
hide solvent
hide ions

# Tampilkan protein sebagai ribbon
show protein cartoons
color protein #cccccc

# Tampilkan ligan sebagai sphere merah
show ligand atoms
style ligand sphere
color ligand #ff0000

# Sembunyikan sidechain protein (opsional, biar bersih)
hide protein atoms

# Pencahayaan soft (opsional)
lighting soft
```

### Troubleshooting Visual

Jika masih ada molekul "nyangkut" yang bukan ligan:

```chimerax
hide ~protein ~ligand
```

---

## Tahap 4: Rekam Video Full HD

### Metode Terbaik (Window Normal, Output HD)

```chimerax
# 1. Pastikan Command Line muncul
ui tool show "Command Line"

# 2. Set ke frame awal
coordset #1 1

# 3. Mulai rekam dengan resolusi output Full HD
movie record size 1920,1080

# 4. Mainkan animasi (ganti 501 dengan jumlah frame-mu)
coordset #1 1,501 pauseFrames 1

# 5. Setelah animasi selesai, encode video
movie encode "E:/Project/herbal-network-pharmacology/videos/264THM_HD.mp4" quality high
```

### Metode Manual (Jika Command Gagal)

1. Ketik: `movie record size 1920,1080` → Enter
2. Klik tombol **Play ▶️** di slider bawah, tunggu sampai habis
3. Ketik: `movie stop` → Enter
4. Ketik: `movie encode "path/video.mp4" quality high` → Enter

---

## Tahap 5: Tips Tambahan

### Resolusi Video

| Resolusi | Command | Untuk |
|----------|---------|-------|
| 720p HD | `size 1280,720` | Presentasi, draft |
| 1080p FHD | `size 1920,1080` | Publikasi jurnal |
| 4K UHD | `size 3840,2160` | Poster, high-end |

### Anti-Aliasing (Jika GPU Kuat)

```chimerax
movie record size 1920,1080 supersample 2
```

> ⚠️ Jangan pakai `supersample 3+` di GPU lemah (MX230, Intel HD) - bisa hang!

### Simpan Session

```chimerax
save session mysession.cxs
```

### Load Session

```chimerax
open mysession.cxs
```

---

## Referensi Path File

### 264THM-PPARG

- Struktur: `kaggle-output-sim/264THM_PPARG/md/md.gro`
- Trajectory: `kaggle-output-sim/264THM_PPARG/md/visual_ready_264THM_fixed.xtc`
- Video Output: `videos/264THM_HD.mp4`

### Luteolin-PDE5A

- Struktur: `kaggle-output-luteolin/Luteolin_PDE5A/md/md.gro`
- Trajectory: `kaggle-output-luteolin/Luteolin_PDE5A/md/visual_ready_Luteolin_fixed.xtc`
- Video Output: `videos/Luteolin_HD.mp4`

---

## Quick Reference Commands

```chimerax
# === SETUP ===
set bg white
hide solvent; hide ions
show protein cartoons; color protein gray
show ligand atoms; style ligand sphere; color ligand red

# === RECORD VIDEO ===
movie record size 1920,1080
coordset #1 1,501 pauseFrames 1
movie encode "output.mp4" quality high

# === UTILITIES ===
view                          # Reset view
turn y 1 360                  # Rotasi 360° sumbu Y
zoom 1.5                      # Zoom in
save snapshot.png supersample 3  # Screenshot HD
```
