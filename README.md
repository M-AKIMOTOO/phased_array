# Phased Array Aligner

山口 2 局データを位相整合して `YAMAGU66` のフェーズドアレイ生データを生成する Rust プログラムです。

## 主な機能
- 幾何遅延 + clock + user 補正をフレーム中心時刻で適用
- 2 局信号を位相整合して合成し、`YAMAGU66_*.raw` を出力
- `YAMAGU66_YAMAGU66_*.cor` を出力（`frinZ` 等で利用可能）
- `ant1/ant2` の振幅スペクトル PNG を出力

## ビルド
```bash
cargo build --release
```

## 実行例
```bash
./target/release/phased_array \
  --schedule r24131a_cor.xml \
  --raw-directory raw \
  --fft 1024 \
  --length 60 \
  --sideband ant1:LSB ant2:LSB
```

## 主なオプション
| オプション | 説明 |
| :--- | :--- |
| `--schedule <XML>` | 相関スケジュール（`.xml`） |
| `--raw-directory <DIR>` | 入力 `.raw` ディレクトリ |
| `--cor-directory <DIR>` | 出力先ディレクトリ |
| `--fft <INT>` | FFT 点数 |
| `--sideband <SB>` | `LSB` / `USB`（`ant1:LSB ant2:USB` 形式可） |
| `--obsfreq <MHZ>` | 基準観測周波数 |
| `--sampling <MHZ>` | サンプリング周波数 |
| `--cpu <INT>` | 使用 CPU スレッド数 |
| `--delay <FLOAT>` | 残留遅延補正（samples） |
| `--rate <FLOAT>` | 残留レート補正（Hz） |

## 出力
`phased_array/YAMAGU66_yyyydddhhmmss/` 配下に以下を出力します。

1. `YAMAGU66_*.raw`（合成後 raw）
2. `YAMAGU66_YAMAGU66_*.cor`（YAMAGU66 の .cor）
3. `YAMAGU66_*_phased_spectrum_amplitude.png`（`phased/ant1/ant2` 振幅スペクトル）
4. `*_geometricdelay.txt`（遅延ログ）

## phased array データの計算式

実装はフレーム中心時刻 `t` で補正量を評価し、`ant1/ant2` を周波数領域で位相整合して合成します。

```text
1) 幾何遅延（baseline = ant2 - ant1）
   tau_g(t) = - ((r2 - r1) · s) / c

2) 相対遅延モデル（clock + user 補正込み）
   clock1(t) = clock1_delay + clock1_rate * t
   clock2(t) = clock2_delay + clock2_rate * t

   tau_rel_no_clock(t) =
       tau_g(t) + coarse_delay + delay/fs + extra_delay_rate*t

   tau_rel(t) = tau_rel_no_clock(t) + (clock2(t) - clock1(t))

   tau1(t) = tau_rel(t) - d_seek
   tau2(t) = 0

3) 整数/小数遅延へ分解
   n1 = round(tau1 * fs),  frac1 = tau1 - n1/fs
   n2 = round(tau2 * fs),  frac2 = tau2 - n2/fs

4) FFT 領域で補正
   X1[k] = FFT( shift(x1[n], n1) )
   X2[k] = FFT( shift(x2[n], n2) )

   X1d[k] = X1[k] * exp(-j*2*pi*f_k*frac1)
   X2d[k] = X2[k] * exp(-j*2*pi*f_k*frac2)

   X1c[k] = X1d[k] * fr_lo1
   X2c[k] = X2d[k] * fr_lo2

5) band-align + 重み付き合成
   Y[k] = w1 * X1c[k] + w2 * X2c_aligned[k]

   y[n] = IFFT(Y[k])  -> 出力ビット数へ量子化 -> YAMAGU66_*.raw

6) YAMAGU66 .cor（1 sector あたり）
   P66[k] = |Y[k]|^2
   Pref   = w1^2 * P1 + w2^2 * P2
   inv    = 1 / (0.5 * Pref * nf * fft^2)
   C66[k] = inv * P66[k]

   C66[k] を実数スペクトル（imag=0）として
   YAMAGU66_YAMAGU66_*.cor に保存
```

`YAMAGU66` の自己相関専用 PNG は出力せず、可視化は `phased/ant1/ant2` を重ねたスペクトル PNG（1枚）を出力します。
