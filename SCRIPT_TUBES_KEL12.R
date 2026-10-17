library(shiny)
library(shinydashboard)
library(leaflet)
library(sf)
library(terra)

# =====================================================
# 1. PATH DATA & PARAMETER
# =====================================================
path_shp <- "D:/WEBGIS/mojokerto_shp.shp"
dir_img  <- "D:/WEBGIS/S2A_MSIL2A_20261001T024141_N0513_R089_T49MFM_20261001T071804.SAFE/S2A_MSIL2A_20261001T024141_N0513_R089_T49MFM_20261001T071804.SAFE/GRANULE/L2A_T49MFM_A058888_20261001T025256/IMG_DATA"
dir_10   <- file.path(dir_img, "R10m")
dir_20   <- file.path(dir_img, "R20m")

file_out  <- "D:/WEBGIS/indeks_scl_mojokerto_v2.tif"      # cache raster NDVI, NDWI, NDBI, SCL
file_stat <- "D:/WEBGIS/statistik_kecamatan_v2.csv"       # cache statistik zonal indeks
file_scl  <- "D:/WEBGIS/scl_kecamatan_v2.csv"             # cache persentase kelas SCL
file_foto <- "D:/WEBGIS/Anggota Kelompok.jpeg"            # foto anggota kelompok
OFFSET    <- 1000   # baseline N0513 (>= 04.00)

# Metode threshold klasifikasi Rendah/Sedang/Tinggi (untuk NDVI, NDWI, NDBI):
#   "kuantil"  = equal weight (persentil 33,3% & 66,7%)
#   "interval" = equal interval (persentil 2-98% dibagi 3)
METODE_THRESHOLD <- "kuantil"

# Rentang tetap legenda indeks (semua indeks)
RG <- c(-1, 1)

# ---- Kelas SCL yang dipakai (Sen2Cor) ----
kode_scl <- c(3, 7, 8, 9, 10)
nama_scl <- c("Bayangan awan (cloud shadow)",
              "Awan probabilitas rendah",
              "Awan probabilitas sedang",
              "Awan probabilitas tinggi",
              "Cirrus")
nama_scl_pendek <- c("Bayangan", "Awan rendah", "Awan sedang", "Awan tinggi", "Cirrus")
warna_scl <- c("#4d4d4d", "#fee391", "#fdae6b", "#e6550d", "#9ecae1")

# ---- Threshold komposit (silakan diubah) ----
TH_AIR    <- 0.0      # Air       : NDWI  >  TH_AIR
TH_VEG    <- 0.3      # Vegetasi  : NDVI  >= TH_VEG
TH_BGN    <- 0.0      # Bangunan  : NDBI  >  TH_BGN
kode_awan <- c(8, 9, 10)   # Awan : SCL awan sedang, awan tinggi, cirrus

kolom_kec <- "WADMKC"
kolom_kab   <- "WADMKK"
kab_dipakai <- c("Mojokerto", "Kabupaten Mojokerto", "Kab. Mojokerto")

JUDUL <- "Analisis Indeks Spektral dan Klasifikasi Tutupan Lahan Kabupaten Mojokerto"

stopifnot(file.exists(path_shp), dir.exists(dir_10), dir.exists(dir_20))
if (!file.exists(file_foto))
  message("Foto tidak ditemukan di: ", file_foto)

# =====================================================
# 2. BATAS ADMINISTRASI (DIFILTER KE MOJOKERTO)
# =====================================================
batas_all <- st_read(path_shp, quiet = TRUE)
if (is.na(st_crs(batas_all))) st_crs(batas_all) <- 4326

cari_kolom <- function(x, nama) {
  hit <- names(x)[toupper(names(x)) == toupper(nama)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

kol_kec <- cari_kolom(batas_all, kolom_kec)
if (is.na(kol_kec))
  stop("Kolom '", kolom_kec, "' tidak ditemukan. Kolom yang ada: ",
       paste(names(batas_all), collapse = ", "))

kol_kab <- cari_kolom(batas_all, kolom_kab)
if (is.na(kol_kab)) {
  message("Kolom '", kolom_kab, "' tidak ada, filter kabupaten dilewati. ",
          "Kolom yang ada: ", paste(names(batas_all), collapse = ", "))
  batas_sf <- batas_all
} else {
  nilai_kab <- trimws(as.character(batas_all[[kol_kab]]))
  message("Isi kolom ", kol_kab, " di shapefile: ",
          paste(sort(unique(nilai_kab)), collapse = " | "))
  pilih <- tolower(nilai_kab) %in% tolower(kab_dipakai)
  if (!any(pilih))
    stop("Tidak ada poligon yang cocok dengan kab_dipakai. Isi kolom ", kol_kab,
         ": ", paste(sort(unique(nilai_kab)), collapse = " | "))
  batas_sf <- batas_all[pilih, ]
}

batas_wgs <- st_transform(batas_sf, 4326)
batas_vec <- vect(batas_wgs)
bb <- st_bbox(batas_wgs)

nama_kec <- as.character(batas_sf[[kol_kec]])
nama_kec[is.na(nama_kec)] <- "Tanpa nama"
n_kec <- length(unique(nama_kec))
message("Jumlah poligon dipakai: ", nrow(batas_sf), " | kecamatan: ", n_kec)

# ---- Luas per kecamatan (hektare), dihitung di UTM 49S ----
luas_poligon_ha <- as.numeric(st_area(st_transform(batas_sf, 32749))) / 10000
luas_kec <- tapply(luas_poligon_ha, nama_kec, sum)

fmt_angka <- function(x, d = 1)
  format(round(x, d), big.mark = ".", decimal.mark = ",", nsmall = d)

# =====================================================
# 3. PROSES BERAT (sekali saja, lalu disimpan ke .tif)
# =====================================================
if (!file.exists(file_out)) {
  
  baca_band <- function(folder, pola) {
    f <- list.files(folder, pattern = pola, full.names = TRUE)
    if (length(f) == 0) stop("Band tidak ditemukan: ", pola)
    rast(f[1])
  }
  
  b03 <- baca_band(dir_10, "_B03_10m\\.jp2$")
  b04 <- baca_band(dir_10, "_B04_10m\\.jp2$")
  b08 <- baca_band(dir_10, "_B08_10m\\.jp2$")
  b11 <- baca_band(dir_20, "_B11_20m\\.jp2$")
  scl <- baca_band(dir_20, "_SCL_20m\\.jp2$")
  
  if (is.na(crs(b04)) || crs(b04) == "") crs(b04) <- "EPSG:32749"
  
  batas_v <- vect(st_transform(batas_sf, crs(b04)))
  if (is.null(intersect(ext(b04), ext(batas_v))))
    stop("Batas tidak overlap dengan tile T49MFM")
  
  cm <- function(r) mask(crop(r, batas_v), batas_v)
  b03 <- cm(b03); b04 <- cm(b04); b08 <- cm(b08)
  b11 <- resample(cm(b11), b04, method = "bilinear")
  scl <- resample(cm(scl), b04, method = "near")      # data kategori: nearest neighbour
  
  ref <- function(r) clamp((r - OFFSET) / 10000, 0, 1)
  B3 <- ref(b03); B4 <- ref(b04); B8 <- ref(b08); B11 <- ref(b11)
  
  idx <- c(
    (B8 - B4)  / (B8 + B4),    # NDVI
    (B3 - B8)  / (B3 + B8),    # NDWI
    (B11 - B8) / (B11 + B8),   # NDBI
    scl                        # SCL (kode kelas)
  )
  names(idx) <- c("NDVI", "NDWI", "NDBI", "SCL")
  writeRaster(idx, file_out, overwrite = TRUE)
}

idx_utm <- rast(file_out)
scl_utm <- idx_utm[["SCL"]]

# =====================================================
# 4. STATISTIK ZONAL PER KECAMATAN (cache CSV)
# =====================================================
nama_indeks <- c("NDVI", "NDWI", "NDBI")   # indeks numerik
urut_indeks <- c("NDVI", "NDBI", "NDWI")   # urutan di tabel analisis
stat_nama   <- c("Mean", "Min", "Max", "Median", "SD")

perlu_hitung <- !file.exists(file_stat) || !file.exists(file_scl)
if (!perlu_hitung) {
  cek <- read.csv(file_stat, stringsAsFactors = FALSE)
  if (!setequal(unique(cek$Kecamatan), unique(nama_kec))) {
    message("Cache statistik tidak cocok dengan shapefile, dihitung ulang.")
    perlu_hitung <- TRUE
  }
}

if (perlu_hitung) {
  message("Menghitung statistik zonal per kecamatan (sekali saja)...")
  
  lv  <- sort(unique(nama_kec))
  bv  <- vect(st_transform(batas_sf, crs(idx_utm)))
  bv$zid <- match(nama_kec, lv)
  zona <- rasterize(bv, idx_utm[[1]], field = "zid")
  kamus <- data.frame(zid = seq_along(lv), Kecamatan = lv)
  
  # ---- statistik NDVI, NDWI, NDBI ----
  fungsi <- list(Mean = "mean", Min = "min", Max = "max",
                 Median = median, SD = sd)
  
  hasil <- do.call(rbind, lapply(names(fungsi), function(f) {
    z <- zonal(idx_utm[[nama_indeks]], zona, fun = fungsi[[f]], na.rm = TRUE)
    do.call(rbind, lapply(nama_indeks, function(nm)
      data.frame(zid = z[[1]], Indeks = nm, Statistik = f, Nilai = z[[nm]])))
  }))
  
  tab <- reshape(hasil, idvar = c("zid", "Indeks"),
                 timevar = "Statistik", direction = "wide")
  names(tab) <- sub("^Nilai\\.", "", names(tab))
  tab   <- merge(tab, kamus, by = "zid")
  tab   <- tab[, c("Kecamatan", "Indeks", "Mean", "Min", "Max", "Median", "SD")]
  tab   <- tab[order(tab$Kecamatan, match(tab$Indeks, nama_indeks)), ]
  write.csv(tab, file_stat, row.names = FALSE)
  
  # ---- persentase piksel tiap kelas SCL per kecamatan ----
  zs <- lapply(kode_scl, function(k)
    zonal(scl_utm == k, zona, fun = "mean", na.rm = TRUE))
  tab_scl <- data.frame(zid = zs[[1]][[1]],
                        sapply(zs, function(z) 100 * z[[2]]))
  names(tab_scl)[-1] <- paste0("SCL_", kode_scl)
  tab_scl <- merge(kamus, tab_scl, by = "zid")
  tab_scl$zid <- NULL
  write.csv(tab_scl, file_scl, row.names = FALSE)
}

stat_tab <- read.csv(file_stat, stringsAsFactors = FALSE)
scl_tab  <- read.csv(file_scl,  stringsAsFactors = FALSE)

ada_data <- tapply(!is.na(stat_tab$Mean), stat_tab$Kecamatan, any)
kec_ada  <- names(ada_data)[ada_data]
stat_tab <- stat_tab[stat_tab$Kecamatan %in% kec_ada, ]
message("Kecamatan yang punya statistik: ", length(kec_ada),
        " dari ", length(unique(nama_kec)))

stat_tab$Luas_ha  <- as.numeric(luas_kec[stat_tab$Kecamatan])
stat_tab$Luas_km2 <- stat_tab$Luas_ha / 100

# ---- Tabel lebar: 1 baris = 1 kecamatan ----
tabel_lebar <- do.call(rbind, lapply(sort(kec_ada), function(k) {
  d   <- stat_tab[stat_tab$Kecamatan == k, ]
  row <- list(Kecamatan = k,
              Luas_ha   = as.numeric(luas_kec[k]),
              Luas_km2  = as.numeric(luas_kec[k]) / 100)
  for (nm in urut_indeks) for (s in stat_nama) {
    v <- d[d$Indeks == nm, s]
    row[[paste(nm, s, sep = "_")]] <- if (length(v)) v[1] else NA_real_
  }
  for (kk in kode_scl) {
    v <- scl_tab[scl_tab$Kecamatan == k, paste0("SCL_", kk)]
    row[[paste0("SCL_", kk)]] <- if (length(v)) v[1] else NA_real_
  }
  as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
}))

fmt <- function(x) formatC(x, format = "f", digits = 3)

warna_idx <- c(NDVI = "#00a65a", NDBI = "#dd4b39", NDWI = "#3c8dbc")

html_tabel_lebar <- function(w) {
  h1 <- paste0(
    "<tr><th rowspan='2'>Kecamatan</th>",
    "<th rowspan='2'>Luas (ha)</th><th rowspan='2'>Luas (km&sup2;)</th>",
    paste0("<th colspan='5' style='background:", warna_idx[urut_indeks], "'>",
           urut_indeks, "</th>", collapse = ""),
    "<th colspan='", length(kode_scl), "' style='background:#605ca8'>SCL (% piksel)</th>",
    "</tr>")
  h2 <- paste0(
    "<tr>",
    paste(rep(paste0("<th class='sub'>", stat_nama, "</th>", collapse = ""),
              length(urut_indeks)), collapse = ""),
    paste0("<th class='sub'>", nama_scl_pendek, "</th>", collapse = ""),
    "</tr>")
  isi <- vapply(seq_len(nrow(w)), function(i) {
    nilai <- unlist(lapply(urut_indeks, function(nm)
      fmt(as.numeric(unlist(w[i, paste(nm, stat_nama, sep = "_")])))))
    nilai_scl <- formatC(as.numeric(unlist(w[i, paste0("SCL_", kode_scl)])),
                         format = "f", digits = 2)
    paste0("<tr><td style='text-align:left'><b>", w$Kecamatan[i], "</b></td>",
           "<td>", fmt_angka(w$Luas_ha[i]), "</td>",
           "<td>", fmt_angka(w$Luas_km2[i], 2), "</td>",
           paste0("<td>", c(nilai, nilai_scl), "</td>", collapse = ""), "</tr>")
  }, character(1))
  HTML(paste0("<div class='tbl-wrap'><table class='tbl-stat'><thead>", h1, h2,
              "</thead><tbody>", paste(isi, collapse = ""), "</tbody></table></div>"))
}

# ---- Popup HTML per kecamatan ----
buat_popup <- function(kec) {
  luas_ha   <- as.numeric(luas_kec[kec])
  teks_luas <- paste0("<div style='margin:6px 0;font-size:15px'>Luas: <b>",
                      fmt_angka(luas_ha), " ha</b> (",
                      fmt_angka(luas_ha / 100, 2), " km&sup2;)</div>")
  
  d <- stat_tab[stat_tab$Kecamatan == kec, ]
  d <- d[match(nama_indeks, d$Indeks), ]
  d <- d[!is.na(d$Indeks), ]
  
  if (nrow(d) == 0)
    return(paste0("<div style='font-size:15px'><b>Kecamatan ", kec, "</b>",
                  teks_luas, "Tidak ada data citra.</div>"))
  
  baris <- paste0(
    "<tr><td><b>", d$Indeks, "</b></td>",
    "<td>", fmt(d$Mean),   "</td>",
    "<td>", fmt(d$Min),    "</td>",
    "<td>", fmt(d$Max),    "</td>",
    "<td>", fmt(d$Median), "</td>",
    "<td>", fmt(d$SD),     "</td></tr>",
    collapse = ""
  )
  
  sc <- scl_tab[scl_tab$Kecamatan == kec, ]
  teks_scl <- if (nrow(sc) > 0) paste0(
    "<div style='margin-top:8px;font-size:14px'><b>Awan dan bayangan (% piksel)</b><br>",
    paste0(nama_scl, ": ",
           fmt_angka(as.numeric(sc[1, paste0("SCL_", kode_scl)]), 2), "%",
           collapse = "<br>"), "</div>") else ""
  
  paste0(
    "<div style='font-size:15px;min-width:460px'><b style='font-size:17px'>Kecamatan ",
    kec, "</b>", teks_luas,
    "<table border='1' cellpadding='9' ",
    "style='border-collapse:collapse;margin-top:6px;text-align:center;width:100%;font-size:14px'>",
    "<tr style='background:#0f766e;color:white'>",
    "<th>Indeks</th><th>Mean</th><th>Min</th><th>Max</th>",
    "<th>Median</th><th>SD</th></tr>", baris, "</table>", teks_scl, "</div>"
  )
}
popup_kec <- vapply(nama_kec, buat_popup, character(1), USE.NAMES = FALSE)

# =====================================================
# 5. SIAPKAN RASTER UNTUK TAMPILAN
# =====================================================
idx_peta <- project(
  aggregate(idx_utm[[nama_indeks]], fact = 4, fun = "mean", na.rm = TRUE),
  "EPSG:4326", method = "bilinear"
)

# SCL: agregasi modus + nearest neighbour (data kategori)
scl_peta <- project(
  aggregate(scl_utm, fact = 4, fun = "modal", na.rm = TRUE),
  "EPSG:4326", method = "near"
)
scl_peta <- resample(scl_peta, idx_peta[[1]], method = "near")

# SCL: hanya 5 kelas yang dipakai (1..5), kelas lain = NA
scl_kat <- classify(scl_peta, cbind(kode_scl, seq_along(kode_scl)), others = NA)
lab_scl <- paste0(nama_scl, " (SCL ", kode_scl, ")")

v_scl <- values(scl_peta, mat = FALSE); v_scl <- v_scl[!is.na(v_scl)]
persen_scl <- vapply(kode_scl, function(k) 100 * mean(v_scl == k), numeric(1))

# ---- Palet kontinu: urutan warna dari NILAI RENDAH ke NILAI TINGGI ----
palet_warna <- list(
  NDVI = c("#a50026", "#f46d43", "#fee08b", "#a6d96a", "#1a9850", "#006837"),
  NDWI = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
  NDBI = c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d")
)
palet_plot <- lapply(palet_warna, function(p) colorRampPalette(p)(100))

# ---- Threshold per indeks (3 kelas: Rendah / Sedang / Tinggi) ----
nama_kelas <- c("Rendah", "Sedang", "Tinggi")
pal_kelas  <- c("#fff7bc", "#fe9929", "#8c2d04")

brk <- lapply(setNames(nama_indeks, nama_indeks), function(nm) {
  v <- values(idx_peta[[nm]], mat = FALSE)
  v <- v[is.finite(v)]
  if (METODE_THRESHOLD == "kuantil") {
    as.numeric(quantile(v, c(0, 1/3, 2/3, 1)))
  } else {
    q <- as.numeric(quantile(v, c(0.02, 0.98)))
    seq(q[1], q[2], length.out = 4)
  }
})

klasifikasi <- function(r, b)
  classify(r, cbind(c(-Inf, b[2], b[3]), c(b[2], b[3], Inf), 1:3),
           include.lowest = TRUE)

kls <- lapply(setNames(nama_indeks, nama_indeks),
              function(nm) klasifikasi(idx_peta[[nm]], brk[[nm]]))

hitung_persen <- function(r, n = 3) {
  v <- values(r, mat = FALSE); v <- v[!is.na(v)]
  100 * as.numeric(table(factor(v, levels = 1:n))) / length(v)
}
fmt_pct <- function(x) paste0(formatC(x, format = "f", digits = 1, decimal.mark = ","), "%")
persen_kelas <- lapply(kls, hitung_persen)

# ---- KOMPOSIT: klasifikasi berbasis threshold ----
# Urutan prioritas: Awan > Air > Vegetasi > Bangunan > Lainnya
#   1 Air, 2 Vegetasi, 3 Bangunan, 4 Awan, 5 Lainnya
r_ndvi <- idx_peta[["NDVI"]]
r_ndwi <- idx_peta[["NDWI"]]
r_ndbi <- idx_peta[["NDBI"]]
r_awan <- scl_peta %in% kode_awan

komp_kls <- ifel(r_awan, 4,
                 ifel(r_ndwi > TH_AIR, 1,
                      ifel(r_ndvi >= TH_VEG, 2,
                           ifel(r_ndbi > TH_BGN, 3, 5))))

nama_komp  <- c("Air", "Vegetasi", "Bangunan", "Awan", "Lainnya")
warna_komp <- c("#3182bd", "#31a354", "#e34a33", "#cccccc", "#fdd49e")
lab_komp <- c(
  paste0("Air (NDWI > ", fmt(TH_AIR), ")"),
  paste0("Vegetasi (NDVI >= ", fmt(TH_VEG), ")"),
  paste0("Bangunan (NDBI > ", fmt(TH_BGN), ")"),
  paste0("Awan (SCL ", paste(kode_awan, collapse = ", "), ")"),
  "Lainnya (tidak memenuhi kriteria)")
persen_komp <- hitung_persen(komp_kls, 5)

teks_metode <- if (METODE_THRESHOLD == "kuantil") {
  paste0("<b>Equal weight (kuantil)</b>: batas kelas diambil dari persentil 33,3% dan 66,7% ",
         "seluruh piksel tiap indeks, sehingga setiap kelas (Rendah, Sedang, Tinggi) ",
         "berisi sekitar sepertiga piksel.")
} else {
  paste0("<b>Equal interval</b>: rentang persentil 2-98% tiap indeks dibagi menjadi 3 selang ",
         "sama lebar.")
}

# ---- Tabel threshold per indeks (dengan persentase piksel) ----
html_threshold <- function() {
  baris <- vapply(urut_indeks, function(nm) {
    b <- brk[[nm]]; pc <- persen_kelas[[nm]]
    paste0("<tr><td><b>", nm, "</b></td>",
           "<td>sampai ", fmt(b[2]), "<br><small>(", fmt_pct(pc[1]), " piksel)</small></td>",
           "<td>", fmt(b[2]), " s.d. ", fmt(b[3]), "<br><small>(", fmt_pct(pc[2]), " piksel)</small></td>",
           "<td>di atas ", fmt(b[3]), "<br><small>(", fmt_pct(pc[3]), " piksel)</small></td></tr>")
  }, character(1))
  HTML(paste0("<div class='tbl-wrap'><table class='tbl-stat'><thead><tr>",
              "<th>Indeks</th><th>Rendah (skor 1)</th><th>Sedang (skor 2)</th><th>Tinggi (skor 3)</th>",
              "</tr></thead><tbody>", paste(baris, collapse = ""),
              "</tbody></table></div>"))
}

# ---- Tabel kelas SCL ----
html_scl <- function() {
  baris <- paste0("<tr><td><b>", kode_scl, "</b></td><td>", nama_scl, "</td><td>",
                  fmt_pct(persen_scl), "</td></tr>", collapse = "")
  HTML(paste0("<div class='tbl-wrap'><table class='tbl-stat'><thead><tr>",
              "<th>Kode SCL</th><th>Kategori</th><th>Persentase piksel</th>",
              "</tr></thead><tbody>", baris, "</tbody></table></div>"))
}

# ---- Tabel threshold komposit ----
html_threshold_komp <- function() {
  aturan <- c(
    paste0("NDWI &gt; ", fmt(TH_AIR)),
    paste0("NDVI &ge; ", fmt(TH_VEG)),
    paste0("NDBI &gt; ", fmt(TH_BGN)),
    paste0("SCL = ", paste(kode_awan, collapse = ", ")),
    "Tidak memenuhi aturan di atas")
  prioritas <- c("3", "4", "5", "1 (tertinggi)", "6 (sisa)")
  baris <- paste0("<tr><td><b>", nama_komp, "</b></td><td>", aturan, "</td><td>",
                  prioritas, "</td><td>", fmt_pct(persen_komp), "</td></tr>", collapse = "")
  HTML(paste0(
    "<div class='tbl-wrap'><table class='tbl-stat'><thead><tr>",
    "<th>Kelas</th><th>Threshold</th><th>Urutan prioritas</th><th>Persentase piksel</th></tr></thead><tbody>",
    baris, "</tbody></table></div>"))
}

# ---- Tabel penentuan kelas komposit: semua 16 kombinasi kondisi ----
html_keputusan <- function() {
  v <- data.frame(
    awan = as.numeric(values(r_awan, mat = FALSE)),
    air  = as.numeric(values(r_ndwi > TH_AIR, mat = FALSE)),
    veg  = as.numeric(values(r_ndvi >= TH_VEG, mat = FALSE)),
    bgn  = as.numeric(values(r_ndbi > TH_BGN, mat = FALSE)))
  v    <- v[complete.cases(v), ]
  kode <- v$awan * 8 + v$air * 4 + v$veg * 2 + v$bgn
  cnt  <- as.numeric(table(factor(kode, levels = 0:15)))
  pct  <- 100 * cnt / sum(cnt)
  
  g <- expand.grid(bgn = 0:1, veg = 0:1, air = 0:1, awan = 0:1)
  g$pct <- pct[g$awan * 8 + g$air * 4 + g$veg * 2 + g$bgn + 1]
  g$kelas <- ifelse(g$awan == 1, 4,
                    ifelse(g$air == 1, 1,
                           ifelse(g$veg == 1, 2,
                                  ifelse(g$bgn == 1, 3, 5))))
  g$alasan <- ifelse(g$awan == 1,
                     "Awan prioritas 1: indeks lain diabaikan",
                     ifelse(g$air == 1,
                            "Bukan awan, NDWI lolos (prioritas 2): Air",
                            ifelse(g$veg == 1,
                                   "Bukan awan & bukan air, NDVI lolos (prioritas 3): Vegetasi",
                                   ifelse(g$bgn == 1,
                                          "Bukan awan, air, vegetasi; NDBI lolos (prioritas 4): Bangunan",
                                          "Tidak ada kriteria yang terpenuhi: Lainnya"))))
  g <- g[order(-g$awan, -g$air, -g$veg, -g$bgn), ]
  
  ya <- function(x) ifelse(x == 1, "<b>Ya</b>", "Tidak")
  teks_putih <- g$kelas %in% c(1, 2, 3)
  baris <- vapply(seq_len(nrow(g)), function(i) {
    paste0("<tr><td>", i, "</td>",
           "<td>", ya(g$awan[i]), "</td><td>", ya(g$air[i]), "</td>",
           "<td>", ya(g$veg[i]),  "</td><td>", ya(g$bgn[i]), "</td>",
           "<td style='background:", warna_komp[g$kelas[i]], ";color:",
           ifelse(teks_putih[i], "white", "black"), ";font-weight:bold'>",
           nama_komp[g$kelas[i]], "</td>",
           "<td style='text-align:left'>", g$alasan[i], "</td>",
           "<td>", fmt_pct(g$pct[i]), "</td></tr>")
  }, character(1))
  
  HTML(paste0(
    "<p>Satu piksel bisa memenuhi lebih dari satu kriteria sekaligus (misalnya NDVI dan NDBI ",
    "sama-sama lolos). Hasil akhirnya ditentukan oleh <b>urutan prioritas</b>: ",
    "Awan, Air, Vegetasi, Bangunan, Lainnya. Tabel berikut memuat semua kombinasi kondisi ",
    "beserta kelas hasilnya dan persentase piksel yang benar-benar jatuh di kombinasi tersebut.</p>",
    "<div class='tbl-wrap'><table class='tbl-stat'><thead><tr>",
    "<th>No</th>",
    "<th>Awan<br><small>(SCL ", paste(kode_awan, collapse = ", "), ")</small></th>",
    "<th>Air<br><small>(NDWI &gt; ", fmt(TH_AIR), ")</small></th>",
    "<th>Vegetasi<br><small>(NDVI &ge; ", fmt(TH_VEG), ")</small></th>",
    "<th>Bangunan<br><small>(NDBI &gt; ", fmt(TH_BGN), ")</small></th>",
    "<th>Kelas hasil</th><th>Alasan penentuan</th><th>Persentase piksel</th>",
    "</tr></thead><tbody>", paste(baris, collapse = ""), "</tbody></table></div>"))
}

# ---- Rumus indeks spektral ----
html_rumus <- function() {
  HTML(paste0(
    "<p><b>Reflektansi</b> (Sentinel-2 L2A): &rho; = (DN &minus; ", OFFSET, ") / 10000, ",
    "dibatasi 0 sampai 1. Band B11 (20 m) di-resample bilinear ke 10 m, ",
    "band SCL (20 m) di-resample nearest neighbour. Nilai tiap indeks berkisar -1 sampai 1.</p>",
    "<div class='tbl-wrap'><table class='tbl-stat'><thead><tr>",
    "<th>Indeks</th><th>Rumus</th><th>Band yang dipakai</th></tr></thead><tbody>",
    "<tr><td><b>NDVI</b></td><td>(B8 &minus; B4) / (B8 + B4)</td>",
    "<td>B8 NIR (10 m), B4 Merah (10 m)</td></tr>",
    "<tr><td><b>NDBI</b></td><td>(B11 &minus; B8) / (B11 + B8)</td>",
    "<td>B11 SWIR-1 (20 m), B8 NIR (10 m)</td></tr>",
    "<tr><td><b>NDWI</b></td><td>(B3 &minus; B8) / (B3 + B8)</td>",
    "<td>B3 Hijau (10 m), B8 NIR (10 m)</td></tr>",
    "<tr><td><b>SCL</b></td><td>Kelas piksel hasil Sen2Cor (bukan hasil hitungan): ",
    "3 = bayangan awan, 7 = awan probabilitas rendah, 8 = awan probabilitas sedang, ",
    "9 = awan probabilitas tinggi, 10 = cirrus</td>",
    "<td>Scene Classification Layer (20 m)</td></tr>",
    "</tbody></table></div>",
    "<p style='margin-top:10px'><b>Statistik zonal per kecamatan</b>: Mean = &Sigma;x / n; ",
    "Min dan Max = nilai terkecil dan terbesar; Median = nilai tengah; ",
    "SD = &radic;[&Sigma;(x &minus; mean)&sup2; / (n &minus; 1)]; luas dari poligon batas (UTM 49S). ",
    "Untuk SCL, yang dihitung adalah persentase piksel tiap kelas terhadap seluruh piksel valid kecamatan.</p>"))
}

# ---- Langkah perhitungan threshold ----
html_langkah <- function() {
  langkah_b <- if (METODE_THRESHOLD == "kuantil") paste0(
    "<li>Hitung dua batas dari persentil: b<sub>1</sub> = P<sub>33,3</sub> dan ",
    "b<sub>2</sub> = P<sub>66,7</sub>. Dengan begitu tiap kelas berisi sekitar 1/3 piksel.</li>") else paste0(
      "<li>Ambil q<sub>2</sub> = P<sub>2</sub> dan q<sub>98</sub> = P<sub>98</sub>; lebar kelas ",
      "w = (q<sub>98</sub> &minus; q<sub>2</sub>) / 3; b<sub>1</sub> = q<sub>2</sub> + w; ",
      "b<sub>2</sub> = q<sub>2</sub> + 2w.</li>")
  HTML(paste0(
    "<h4>A. Klasifikasi Rendah / Sedang / Tinggi (NDVI, NDWI, NDBI)</h4>",
    "<p>Metode yang dipakai: ", teks_metode, "</p><ol>",
    "<li>Kumpulkan semua nilai piksel (x) tiap indeks pada raster tampilan 40 m, tanpa NA.</li>",
    langkah_b,
    "<li>Klasifikasi: <b>Rendah</b> jika x &le; b<sub>1</sub>; ",
    "<b>Sedang</b> jika b<sub>1</sub> &lt; x &le; b<sub>2</sub>; ",
    "<b>Tinggi</b> jika x &gt; b<sub>2</sub>.</li></ol>",
    "<h4>B. Komposit (Air, Vegetasi, Bangunan, Awan)</h4>",
    "<p>Tiap piksel dicek berurutan dengan threshold tetap; aturan pertama yang terpenuhi dipakai:</p><ol>",
    "<li><b>Awan</b> jika SCL = ", paste(kode_awan, collapse = ", "), ".</li>",
    "<li><b>Air</b> jika NDWI &gt; ", fmt(TH_AIR), ".</li>",
    "<li><b>Vegetasi</b> jika NDVI &ge; ", fmt(TH_VEG), ".</li>",
    "<li><b>Bangunan</b> jika NDBI &gt; ", fmt(TH_BGN), ".</li>",
    "<li>Selain itu masuk <b>Lainnya</b> (misalnya lahan terbuka atau vegetasi jarang).</li>",
    "</ol>"))
}

panah_utara_html <- paste0(
  "<div style='background:rgba(255,255,255,.88);padding:4px 9px;border-radius:6px;",
  "border:1px solid #888;text-align:center'>",
  "<div style='font-weight:bold;font-size:16px;line-height:16px'>N</div>",
  "<svg width='26' height='40' viewBox='0 0 26 40'>",
  "<polygon points='13,2 24,38 13,29 2,38' fill='#222'/></svg></div>")

# ---- Arah mata angin (utara) untuk plot dasar R ----
tambah_arah_utara <- function() {
  u <- par("usr"); w <- u[2] - u[1]; h <- u[4] - u[3]
  xn <- u[2] - 0.07 * w; yt <- u[4] - 0.06 * h; yb <- u[4] - 0.17 * h
  rect(xn - 0.04 * w, yb - 0.02 * h, xn + 0.04 * w, yt + 0.06 * h,
       col = adjustcolor("white", 0.75), border = NA, xpd = NA)
  polygon(c(xn, xn - 0.015 * w, xn, xn + 0.015 * w),
          c(yt, yb, yb + 0.035 * h, yb), col = "black", border = "black", xpd = NA)
  text(xn, yt + 0.03 * h, "N", font = 2, cex = 1.2, xpd = NA)
}

# ---- Skala batang (km) untuk plot dasar R; raster berkoordinat lon/lat ----
tambah_skala <- function(r) {
  e   <- ext(r)
  lat <- (ymin(e) + ymax(e)) / 2
  lebar_km <- (xmax(e) - xmin(e)) * 111.32 * cos(lat * pi / 180)
  pilihan  <- c(1, 2, 5, 10, 20, 25, 50)
  d <- pilihan[which.min(abs(pilihan - lebar_km / 4))]
  sbar(d, xy = "bottomleft", type = "bar", divs = 2, below = "km",
       labels = c(0, d / 2, d), cex = 1.1)
}

# ---- Colorbar -1 s.d. +1 (di margin kanan plot) ----
tambah_colorbar <- function(rg, pal, judul) {
  op <- par(xpd = NA); on.exit(par(op))
  u <- par("usr"); w <- u[2] - u[1]; h <- u[4] - u[3]
  x0 <- u[2] + 0.03 * w; x1 <- x0 + 0.05 * w
  y0 <- u[3] + 0.12 * h; y1 <- u[4] - 0.12 * h
  n  <- length(pal); ys <- seq(y0, y1, length.out = n + 1)
  rect(x0, ys[-(n + 1)], x1, ys[-1], col = pal, border = NA)
  rect(x0, y0, x1, y1, border = "black", lwd = 1.5)
  tk <- seq(rg[1], rg[2], by = 0.5)
  yt <- y0 + (tk - rg[1]) / (rg[2] - rg[1]) * (y1 - y0)
  segments(x1, yt, x1 + 0.012 * w, yt)
  text(x1 + 0.018 * w, yt, sprintf("%+.1f", tk), adj = 0, cex = 1.2)
  text((x0 + x1) / 2, y1 + 0.05 * h, judul, font = 2, cex = 1.3)
}

# =====================================================
# 6. UI
# =====================================================
css_tema <- "
/* ---------- Header ---------- */
.skin-blue .main-header .logo{background:#0f766e;color:#fff;font-weight:700;}
.skin-blue .main-header .logo:hover{background:#115e59;}
.skin-blue .main-header .navbar{background:#0f766e;}
.skin-blue .main-header .navbar .sidebar-toggle:hover{background:#115e59;}
.main-header .navbar-custom-menu{float:left !important;}
.judul-header{display:block;color:#fff;font-size:20px;font-weight:600;line-height:50px;
  padding:0 15px;white-space:nowrap;}
@media (max-width:991px){
  .judul-header{font-size:13px;line-height:18px;padding-top:6px;white-space:normal;}
}
/* ---------- Sidebar ---------- */
.skin-blue .wrapper,.skin-blue .main-sidebar,.skin-blue .left-side{background-color:#134e4a;}
.skin-blue .sidebar a{color:#ccfbf1;}
.skin-blue .sidebar-menu>li>a{border-left:3px solid transparent;}
.skin-blue .sidebar-menu>li:hover>a,.skin-blue .sidebar-menu>li.active>a,
.skin-blue .sidebar-menu>li.menu-open>a{color:#fff;background:#115e59;border-left-color:#f59e0b;}
.skin-blue .sidebar-menu>li>.treeview-menu{margin:0;background:#0f3d3a;}
.skin-blue .treeview-menu>li>a{color:#99f6e4;}
.skin-blue .treeview-menu>li.active>a,.skin-blue .treeview-menu>li>a:hover{color:#fbbf24;}
/* ---------- Konten & box ---------- */
.content-wrapper,.right-side{background:#f0fdfa;}
.box.box-solid.box-primary{border:1px solid #0f766e;}
.box.box-solid.box-primary>.box-header{background:#0f766e;background-color:#0f766e;color:#fff;}
.box.box-solid.box-info{border:1px solid #0f766e;}
.box.box-solid.box-info>.box-header{background:#b45309;background-color:#b45309;color:#fff;}
.btn-default{border-color:#0f766e;color:#0f766e;}
.btn-default:hover{background:#0f766e;color:#fff;}
/* ---------- Tabel ---------- */
.tbl-wrap{overflow-x:auto;}
.tbl-stat{border-collapse:collapse;width:100%;font-size:13px;text-align:center;}
.tbl-stat th,.tbl-stat td{border:1px solid #999;padding:6px 9px;white-space:nowrap;}
.tbl-stat th{background:#0f766e;color:#fff;}
.tbl-stat th.sub{background:#ccfbf1;color:#222;}
.tbl-stat tbody tr:nth-child(even){background:#f0fdfa;}
/* ---------- Peta ---------- */
.leaflet-popup-content{font-size:15px;}
.leaflet-control.legend, .info.legend{font-size:15px;line-height:22px;padding:8px 12px;}
.info.legend i{width:22px;height:22px;}
.info.legend .legend-title, .info.legend div{font-size:15px;}
/* ---------- Identitas ---------- */
.identitas h2{margin-top:0;color:#0f766e;font-weight:800;letter-spacing:.5px;}
.identitas h3{margin-top:0;color:#0f766e;font-weight:700;}
.identitas p{font-size:18px;}
#foto_anggota img{max-width:100%;height:auto;border-radius:10px;border:3px solid #0f766e;}
/* ---------- Sumber data ---------- */
.sumber-data li{margin-bottom:6px;font-size:14px;}
"

ui <- dashboardPage(
  dashboardHeader(
    title = "WebGIS", titleWidth = 300,
    tags$li(class = "dropdown", tags$span(class = "judul-header", JUDUL))
  ),
  
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "menu",
      menuItem("Identitas", tabName = "m_id", icon = icon("users"), selected = TRUE),
      menuItem(
        text = "Maps", tabName = "maps", icon = icon("globe"),
        menuSubItem("Peta Interaktif",     tabName = "m_osm",  icon = icon("map")),
        menuSubItem("Hasil Plot Indeks",   tabName = "m_plot", icon = icon("chart-area")),
        menuSubItem("Hasil Komposit",      tabName = "m_komp", icon = icon("layer-group"))
      ),
      menuItem("Analisis", tabName = "analisis", icon = icon("table"),
               menuSubItem("Tabel Statistik Kecamatan", tabName = "m_stat", icon = icon("table")),
               menuSubItem("Threshold dan Perhitungan", tabName = "m_thr",  icon = icon("calculator")))
    )
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML(css_tema))),
    tabItems(
      
      # ---- Identitas ----
      tabItem(
        tabName = "m_id",
        box(
          width = 12, title = "Identitas", status = "primary", solidHeader = TRUE,
          div(class = "identitas",
              h2("KELOMPOK 12 PEMROGRAMAN SPASIAL KELAS A"),
              imageOutput("foto_anggota", height = "auto"),
              br(),
              h3("Lokasi Studi"),
              p(paste0("Kabupaten Mojokerto, Jawa Timur dengan ", n_kec, " Kecamatan")),
              leafletOutput("map_id", height = "480px"),
              br(),
              h3("Alasan Pemilihan Lokasi"),
              p(paste0("Kabupaten Mojokerto yang dipilih karena memiliki karakteristik tutupan lahan ",
                       "yang beragam, seperti kawasan permukiman, industri, tambak, dan lahan pertanian. ",
                       "Keberagaman tersebut memungkinkan identifikasi perbedaan karakteristik spektral ",
                       "objek menggunakan citra Sentinel-2."))
          )
        )
      ),
      
      # ---- Peta interaktif ----
      tabItem(
        tabName = "m_osm",
        fluidRow(
          column(
            width = 3,
            box(
              width = 12, title = "Pengaturan", status = "primary", solidHeader = TRUE,
              radioButtons(
                "indeks", "Pilih layer:",
                choices = c(
                  "NDVI" = "NDVI", "NDWI" = "NDWI", "NDBI" = "NDBI",
                  "SCL (awan dan bayangan)" = "SCL",
                  "Komposit berdasarkan threshold" = "komp"
                )
              ),
              sliderInput("op", "Opacity raster", 0, 1, 0.8, 0.1),
              checkboxInput("tampil_nama", "Tampilkan nama kecamatan di peta", TRUE),
              helpText("Klik poligon kecamatan untuk melihat luas dan tabel statistik zonal. ",
                       "Batas administrasi bisa dinyalakan atau dimatikan lewat kontrol layer di peta."),
              hr(),
              downloadButton("dl", "Download GeoTIFF")
            ),
            box(
              width = 12, title = "Sumber Data", status = "info", solidHeader = TRUE,
              tags$ul(class = "sumber-data",
                      tags$li(tags$b("Batas administrasi (shapefile): "),
                              "BIG Inageoportal (Badan Informasi Geospasial)."),
                      tags$li(tags$b("Citra satelit: "),
                              "Sentinel-2 Level-2A dari Copernicus (ESA).")
              )
            )
          ),
          column(width = 9, leafletOutput("map", height = "85vh"))
        )
      ),
      
      # ---- Hasil plot 4 panel ----
      tabItem(
        tabName = "m_plot",
        box(
          width = 12, title = "Hasil Plot NDVI, NDWI, NDBI, dan SCL",
          status = "primary", solidHeader = TRUE,
          plotOutput("plot4", height = "850px")
        )
      ),
      
      # ---- Hasil komposit ----
      tabItem(
        tabName = "m_komp",
        box(width = 12, title = "Analisis Tutupan Lahan di Kabupaten Mojokerto",
            status = "primary", solidHeader = TRUE,
            plotOutput("plot_komp", height = "750px"))
      ),
      
      # ---- Analisis: tabel statistik + peringkat + diagram ----
      tabItem(
        tabName = "m_stat",
        box(
          width = 12, title = "Statistik Zonal per Kecamatan",
          status = "primary", solidHeader = TRUE,
          fluidRow(
            column(4, selectInput("pilih_kec", "Kecamatan:",
                                  c("Semua", sort(kec_ada)))),
            column(4, br(), downloadButton("dl_stat", "Download CSV"))
          ),
          uiOutput("tabel_stat"),
          helpText("Satu baris = satu kecamatan; tiap indeks memuat Mean, Min, Max, Median, dan SD.",
                   "Kolom SCL berisi persentase piksel tiap kelas terhadap seluruh piksel valid kecamatan.",
                   "Statistik dihitung dari raster 10 m (bukan versi untuk tampilan peta).",
                   "Luas dihitung dari poligon shapefile (UTM 49S).")
        ),
        box(
          width = 12, title = "Peringkat Kecamatan dan Diagram Statistik",
          status = "primary", solidHeader = TRUE,
          radioButtons("idx_diag", "Pilih indeks:",
                       choices = c("NDVI", "NDBI", "NDWI"), inline = TRUE),
          uiOutput("ringkas_idx"),
          hr(),
          plotOutput("diagram_stat", height = "auto"),
          helpText("Peringkat berdasarkan nilai rata-rata (Mean) indeks per kecamatan. ",
                   "Pada diagram, kecamatan diurutkan dari rata-rata tertinggi (paling atas) ",
                   "ke terendah (paling bawah). Garis hitam menandai nilai 0.")
        )
      ),
      
      # ---- Analisis: threshold dan perhitungan ----
      tabItem(
        tabName = "m_thr",
        box(width = 12, title = "1. Rumus Indeks Spektral dan SCL",
            status = "primary", solidHeader = TRUE, html_rumus()),
        box(width = 12, title = "2. Cara Menghitung Threshold dan Komposit",
            status = "primary", solidHeader = TRUE, html_langkah()),
        box(width = 12, title = "3. Tabel Threshold Kategori per Indeks",
            status = "primary", solidHeader = TRUE, html_threshold(),
            helpText("Persentase piksel dihitung dari raster tampilan 40 m seluruh Kabupaten Mojokerto.")),
        box(width = 12, title = "4. Tabel Kelas SCL",
            status = "primary", solidHeader = TRUE, html_scl()),
        box(width = 12, title = "5. Tabel Threshold Komposit",
            status = "primary", solidHeader = TRUE, html_threshold_komp()),
        box(width = 12, title = "6. Tabel Penentuan Kelas Komposit (Semua Kombinasi Kondisi)",
            status = "primary", solidHeader = TRUE, html_keputusan())
      )
    )
  )
)

# =====================================================
# 7. SERVER
# =====================================================
server <- function(input, output, session) {
  
  # ---- Foto anggota kelompok ----
  output$foto_anggota <- renderImage({
    if (!file.exists(file_foto)) return(NULL)
    list(src = file_foto, contentType = "image/jpeg", alt = "Foto anggota kelompok")
  }, deleteFile = FALSE)
  
  output$map <- renderLeaflet({
    leaflet() %>%
      addMapPane("raster", zIndex = 450) %>%
      addMapPane("batas",  zIndex = 460) %>%
      addProviderTiles("OpenStreetMap", group = "OSM") %>%
      addProviderTiles("Esri.WorldImagery", group = "Satelit") %>%
      addLayersControl(
        baseGroups    = c("OSM", "Satelit"),
        overlayGroups = c("Indeks", "Batas"),
        options       = layersControlOptions(collapsed = FALSE)
      ) %>%
      addScaleBar(position = "bottomleft",
                  options = scaleBarOptions(metric = TRUE, imperial = FALSE)) %>%
      addControl(html = panah_utara_html, position = "topleft") %>%
      fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
  })
  
  observe({
    jenis <- input$indeks
    
    proxy <- leafletProxy("map") %>%
      clearGroup("Indeks") %>%
      clearGroup("Batas") %>%
      clearPopups() %>%
      removeControl("legenda")
    
    if (jenis == "komp") {
      pal_k <- colorFactor(warna_komp, levels = 1:5, na.color = "transparent")
      proxy %>%
        addRasterImage(komp_kls, colors = pal_k, opacity = input$op,
                       project = FALSE, group = "Indeks",
                       options = gridOptions(pane = "raster")) %>%
        addLegend("bottomright", colors = warna_komp, labels = lab_komp,
                  title = "Komposit threshold", opacity = 1, layerId = "legenda")
      
    } else if (jenis == "SCL") {
      pal_s <- colorFactor(warna_scl, levels = 1:5, na.color = "transparent")
      proxy %>%
        addRasterImage(scl_kat, colors = pal_s, opacity = input$op,
                       project = FALSE, group = "Indeks",
                       options = gridOptions(pane = "raster")) %>%
        addLegend("bottomright", colors = warna_scl, labels = lab_scl,
                  title = "SCL", opacity = 1, layerId = "legenda")
      
    } else {
      r   <- idx_peta[[jenis]]
      pal <- colorNumeric(palet_warna[[jenis]], domain = RG, na.color = "transparent")
      proxy %>%
        addRasterImage(r, colors = pal, opacity = input$op, project = FALSE,
                       group = "Indeks", options = gridOptions(pane = "raster")) %>%
        addLegend("bottomright", pal = pal, values = RG, title = jenis,
                  labFormat = function(type, cuts, p) sprintf("%+.1f", cuts),
                  opacity = 1, layerId = "legenda")
    }
    
    # ---- Batas administrasi (selalu dipasang; hidup/mati lewat kontrol layer) ----
    opsi_label <- if (isTRUE(input$tampil_nama)) {
      gaya <- list(
        "font-size"   = "12px",
        "font-weight" = "bold",
        "color"       = "#111111",
        "text-shadow" = "-1px -1px 2px #fff, 1px -1px 2px #fff, -1px 1px 2px #fff, 1px 1px 2px #fff"
      )
      labelOptions(noHide = TRUE, direction = "center",
                   textOnly = TRUE, style = gaya)
    } else {
      labelOptions(direction = "auto", sticky = TRUE)
    }
    
    leafletProxy("map") %>%
      addPolygons(
        data = batas_wgs, group = "Batas",
        layerId = as.character(seq_len(nrow(batas_wgs))),
        color = "black", weight = 2,
        fillColor = "#000000", fillOpacity = 0.01,
        label = nama_kec,
        labelOptions = opsi_label,
        highlightOptions = highlightOptions(weight = 4, color = "yellow",
                                            bringToFront = TRUE),
        options = pathOptions(pane = "batas")
      )
  })
  
  observeEvent(input$map_shape_click, {
    klik <- input$map_shape_click
    if (is.null(klik$id)) return()
    i <- suppressWarnings(as.integer(klik$id))
    if (is.na(i) || i < 1 || i > length(popup_kec)) return()
    
    leafletProxy("map") %>%
      clearPopups() %>%
      addPopups(lng = klik$lng, lat = klik$lat, popup = popup_kec[i],
                options = popupOptions(maxWidth = 700, minWidth = 480))
  })
  
  output$dl <- downloadHandler(
    filename = function() {
      nm <- if (input$indeks %in% c(nama_indeks, "SCL")) input$indeks else "komposit_threshold"
      paste0(nm, "_mojokerto.tif")
    },
    content = function(file) {
      if (input$indeks %in% c(nama_indeks, "SCL"))
        writeRaster(idx_utm[[input$indeks]], file, overwrite = TRUE)
      else
        writeRaster(komp_kls, file, overwrite = TRUE)
    }
  )
  
  # ---- Plot 2x2: NDVI, NDWI, NDBI (skala -1 s.d. 1) + SCL ----
  output$plot4 <- renderPlot({
    par(mfrow = c(2, 2))
    for (nm in nama_indeks) {
      plot(idx_peta[[nm]], col = palet_plot[[nm]], range = RG,
           legend = FALSE, axes = TRUE, main = nm,
           mar = c(3, 3, 3, 8))
      plot(batas_vec, add = TRUE, border = "black", lwd = 0.8)
      tambah_arah_utara()
      tambah_colorbar(RG, palet_plot[[nm]], nm)
    }
    plot(scl_kat, breaks = seq(0.5, 5.5, by = 1), col = warna_scl,
         legend = FALSE, axes = TRUE, main = "SCL",
         mar = c(3, 3, 3, 8))
    plot(batas_vec, add = TRUE, border = "black", lwd = 0.8)
    tambah_arah_utara()
    legend("bottomright", title = "SCL", legend = lab_scl,
           fill = warna_scl, bg = "white", cex = 0.8)
  })
  
  # ---- Komposit: judul, legenda, arah utara, dan skala ----
  output$plot_komp <- renderPlot({
    par(mfrow = c(1, 1))
    plot(komp_kls, breaks = seq(0.5, 5.5, by = 1), col = warna_komp,
         legend = FALSE, axes = TRUE,
         main = "Analisis Tutupan Lahan di Kabupaten Mojokerto",
         cex.main = 1.6, mar = c(3, 3, 4, 3))
    plot(batas_vec, add = TRUE, border = "black", lwd = 0.8)
    tambah_arah_utara()
    tambah_skala(komp_kls)
    legend("bottomright", title = "Komposit threshold", legend = lab_komp,
           fill = warna_komp, bg = "white", cex = 1.1)
  })
  
  tab_tampil <- reactive({
    d <- tabel_lebar
    if (input$pilih_kec != "Semua") d <- d[d$Kecamatan == input$pilih_kec, ]
    d
  })
  
  output$tabel_stat <- renderUI(html_tabel_lebar(tab_tampil()))
  
  output$dl_stat <- downloadHandler(
    filename = function() "statistik_kecamatan_mojokerto.csv",
    content  = function(file) write.csv(tab_tampil(), file, row.names = FALSE)
  )
  
  # ---- Kecamatan tertinggi dan terendah (berdasarkan rata-rata indeks) ----
  output$ringkas_idx <- renderUI({
    nm <- input$idx_diag
    m  <- tabel_lebar[[paste0(nm, "_Mean")]]
    k  <- tabel_lebar$Kecamatan
    hi <- which.max(m)
    lo <- which.min(m)
    if (length(hi) == 0 || length(lo) == 0) return(NULL)
    
    kotak <- function(judul, kec, nilai, warna) paste0(
      "<div style='flex:1;min-width:260px;background:", warna,
      ";color:#fff;padding:14px 18px;border-radius:8px'>",
      "<div style='font-size:14px;opacity:.9'>", judul, "</div>",
      "<div style='font-size:24px;font-weight:700'>", kec, "</div>",
      "<div style='font-size:15px'>Rata-rata ", nm, " = ", fmt(nilai), "</div></div>")
    
    HTML(paste0(
      "<div style='display:flex;gap:14px;flex-wrap:wrap'>",
      kotak(paste("Kecamatan dengan", nm, "tertinggi"), k[hi], m[hi], "#0f766e"),
      kotak(paste("Kecamatan dengan", nm, "terendah"),  k[lo], m[lo], "#b45309"),
      "</div>"))
  })
  
  # ---- Diagram Max, Min, Mean, SD, Median per kecamatan untuk indeks terpilih ----
  output$diagram_stat <- renderPlot({
    nm <- input$idx_diag
    d  <- tabel_lebar
    d  <- d[order(d[[paste0(nm, "_Mean")]], na.last = FALSE), ]   # tertinggi di paling atas
    
    stat_plot  <- c("Max", "Min", "Mean", "SD", "Median")
    warna_stat <- c(Max = "#d73027", Min = "#4575b4", Mean = "#1a9850",
                    SD = "#636363", Median = "#7b3294")
    
    par(mfrow = c(1, 5), oma = c(0, 0, 2.5, 0))
    for (i in seq_along(stat_plot)) {
      s <- stat_plot[i]
      v <- d[[paste(nm, s, sep = "_")]]
      v[is.na(v)] <- 0
      par(mar = c(4, if (i == 1) 10 else 1, 3, 1))
      
      xl <- range(c(v, 0))
      if (diff(xl) == 0) xl <- xl + c(-1, 1)
      xl <- xl + c(-1, 1) * 0.25 * diff(xl)
      
      bp <- barplot(v, names.arg = if (i == 1) d$Kecamatan else NULL,
                    horiz = TRUE, las = 1, col = warna_stat[[s]], xlim = xl,
                    cex.names = 0.85, main = s, xlab = paste(nm, s))
      text(v, bp, sprintf("%.2f", v), pos = ifelse(v >= 0, 4, 2),
           cex = 0.75, xpd = NA)
      abline(v = 0, lwd = 1.5)
    }
    mtext(paste("Statistik", nm, "per kecamatan"), outer = TRUE, font = 2, cex = 1.2)
  }, height = function() max(450, 26 * nrow(tabel_lebar) + 150))
  
  # ---- Peta shapefile di halaman Identitas ----
  output$map_id <- renderLeaflet({
    pal_id <- colorFactor("Set3", domain = nama_kec)
    leaflet(batas_wgs) %>%
      addProviderTiles("OpenStreetMap") %>%
      addPolygons(color = "black", weight = 1.5,
                  fillColor = ~pal_id(nama_kec), fillOpacity = 0.55,
                  label = nama_kec,
                  labelOptions = labelOptions(noHide = TRUE, direction = "center",
                                              textOnly = TRUE,
                                              style = list("font-size" = "11px",
                                                           "font-weight" = "bold")),
                  highlightOptions = highlightOptions(weight = 3, color = "yellow",
                                                      bringToFront = TRUE)) %>%
      addScaleBar(position = "bottomleft",
                  options = scaleBarOptions(metric = TRUE, imperial = FALSE)) %>%
      addControl(html = panah_utara_html, position = "topleft") %>%
      fitBounds(bb[["xmin"]], bb[["ymin"]], bb[["xmax"]], bb[["ymax"]])
  })
}

# =====================================================
# 8. JALANKAN
# =====================================================
runApp(shinyApp(ui, server), launch.browser = TRUE)