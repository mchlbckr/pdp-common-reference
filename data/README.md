# Data

No raw data is redistributed in this repository. Every dataset is obtained from
its original source with the commands below and verified against the SHA-256
checksums given here. All paths are relative to the repository root, and the
pipeline expects exactly these locations.

Run everything from the repository root, then verify:

```sh
shasum -a 256 -c data/checksums.sha256
```

## Wine Quality (red subset)

The red-wine subset of the UCI Wine Quality data (Cortez et al., 2009). Used
only for the held-out predictive demonstration in S16.

* Source: <https://archive.ics.uci.edu/dataset/186/wine+quality>
* Licence: CC BY 4.0
* Rows: 1,599 observations plus a header; 11 physicochemical predictors and the
  recorded sensory-quality score.
* Target path: `data/winequality-red.csv`
* SHA-256: `4a402cf041b025d4566d954c3b9ba8635a3a8a01e039005d97d6a710278cf05e`

```sh
mkdir -p data
curl -L -o /tmp/wine-quality.zip \
  https://archive.ics.uci.edu/static/public/186/wine+quality.zip
unzip -o -j /tmp/wine-quality.zip winequality-red.csv -d data
```

## Bank Marketing

Held-out mixed-type classification application. Also supplies the controlled
off-support stress test reported in the manuscript.

* Source: <https://archive.ics.uci.edu/dataset/222/bank+marketing>
* DOI: `10.24432/C5K306`
* Licence: CC BY 4.0
* Rows: 45,211 records.
* Target path: `data/uci/bank-extract/bank-full.csv`
* SHA-256: `d1513ec63b385506f7cfce9f2c5caa9fe99e7ba4e8c3fa264b3aaf0f849ed32d`

The UCI archive nests a second archive: `bank+marketing.zip` contains
`bank.zip`, which in turn holds `bank-full.csv`.

```sh
mkdir -p data/uci/bank-extract
curl -L -o data/uci/bank-marketing.zip \
  https://archive.ics.uci.edu/static/public/222/bank+marketing.zip
unzip -o data/uci/bank-marketing.zip -d data/uci/bank-extract
unzip -o data/uci/bank-extract/bank.zip -d data/uci/bank-extract
```

## Online Shoppers Purchasing Intention

Held-out mixed-type classification application used as a negative control.

* Source: <https://archive.ics.uci.edu/dataset/468/online+shoppers+purchasing+intention+dataset>
* DOI: `10.24432/C5F88Q`
* Licence: CC BY 4.0
* Rows: 12,330 sessions.
* Target path: `data/uci/online_shoppers_intention.csv`
* SHA-256: `b3055ee355f59134d851d32641183cb4a8b45def7124d2f50442a042f358e0d9`

```sh
mkdir -p data/uci
curl -L -o data/uci/online-shoppers.zip \
  https://archive.ics.uci.edu/static/public/468/online+shoppers+purchasing+intention+dataset.zip
unzip -o data/uci/online-shoppers.zip -d data/uci
```

## California Housing

The 20,640 California census-block observations introduced by Pace and Barry
(1997) and distributed through StatLib. This is the confirmatory application of
the manuscript, analysed over 30 independently seeded splits in
`simulations/s33_california_longitude_application.R`.

* Source: <http://lib.stat.cmu.edu/datasets/houses.zip>
* Original article DOI: `10.1016/S0167-7152(96)00140-X`
* Rows: 20,640; eight derived or original predictors and median house value.
* Target path: `data/pilot/california/CaliforniaHousing/cal_housing.data`
* SHA-256 of `cal_housing.data`: `a2b44dcc7e1c9091b43c4f458b25cbd1f7e8068d08edd26bc4703108d3fa6214`
* SHA-256 of `houses.zip`: `aaa5c9a6afe2225cc2aed2723682ae403280c4a3695a2ddda4ffb5d8215ea681`

The directory is named `pilot/` for historical reasons. The path is hard-coded
in `exploration/s33_real_data_screen.R` and must be kept as given.

```sh
mkdir -p data/pilot/california
curl -L -o data/pilot/california/houses.zip \
  http://lib.stat.cmu.edu/datasets/houses.zip
unzip -o data/pilot/california/houses.zip -d data/pilot/california
```

StatLib is occasionally unavailable. The same data ships with `scikit-learn` as
`sklearn.datasets.fetch_california_housing`, but that version is preprocessed
and is **not** byte-identical to `cal_housing.data`; use it only as a fallback
for inspection, never for reproducing the reported numbers.
