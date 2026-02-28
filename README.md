# Phased Array Aligner

山口 32m-34m 干渉計のデータを相関処理し、フェーズドアーレイ（結合型干渉計）出力を生成するための Rust プログラムです。

## 主な機能
- **高精度幾何学的補正**: 地球自転に伴う幾何学的遅延を 1秒ごとに更新し、フレーム単位で精密に補正（Delay/Rate/Accel）します。
- **フリンジ解析 (`--fringe`)**: 相互相関（XCF）を計算し、遅延ピーク検出やスペクトル・位相のプロットを生成します。
- **フェーズドアーレイ合成**: 2つのアンテナ信号の位相を完全に同期させて加算し、1つの強力な単一アンテナデータ（YAMAGU66）として出力します。
- **柔軟な個別設定**: アンテナごとに異なるサイドバンド（LSB/USB）、ビット数、量子化レベル、ビットシャッフルを指定可能です。
- **自動スペクトル反転**: LSB データをデコード時に自動反転し、内部計算と出力をすべて周波数増加方向（USB）に統一します。

## ビルド方法
Rust の開発環境（Cargo）が必要です。常に `--release` フラグを付けてコンパイルしてください。
```bash
cargo build --release
```

## 使い方

### 基本的な実行例
`.ifile`（または XML スケジュール）とデータディレクトリを指定して実行します。
```bash
./target/release/phased_array --ifile 3C345.ifile --data raw --fft 1024 --length 60 --fringe --sideband ant1:LSB ant2:LSB
```

### 主要なオプション

| オプション | 説明 |
| :--- | :--- |
| `--ifile <PATH>` | 天体座標やアンテナ位置が記された設定ファイル（.ifile または .xml）。 |
| `--data <DIR>` | `.raw` データファイルが格納されているディレクトリ。 |
| `--fft <INT>` | FFT ポイント数。相関探索には 1024、高分解能合成には 8192 以上を推奨。 |
| `--fringe` | フリンジ解析を実行し、`results` ディレクトリに各種プロットを出力。 |
| `--sideband <SB>` | `LSB` または `USB`。 `ant1:LSB ant2:USB` のように個別指定も可能。 |
| `--vsrec` | `--level` と `--shuffle` を VSREC 固定値へ強制上書き（2-bit 専用）。 |
| `--obsfreq <MHZ>` | 観測帯域の下限周波数（MHz）。ドップラー補正の基準（0 Hz）として使用。 |
| `--sampling <MHZ>` | サンプリング周波数（MHz、デフォルト: 1024.0）。 |
| `--cpu <INT>` | 使用する CPU コア数。 |
| `--delay <FLOAT>` | 残留遅延の追加補正（サンプル単位）。 |
| `--rate <FLOAT>` | 残留レートの追加補正（Hz単位）。 |

## 出力ファイル

実行後、`phased_array/YAMAGU66_yyyydddhhmmss/` ディレクトリに以下のファイルが生成されます。

1.  **`.raw` ファイル**: フェーズドアーレイ合成後の 2-bit (or 4-bit) 生データ。
2.  **`.cor` ファイル**: `frinZ` 等で解析可能な相関スペクトルデータ（4種類）。
    - `YAMAGU32_YAMAGU32_...`: アンテナ1 自己相関
    - `YAMAGU34_YAMAGU34_...`: アンテナ2 自己相関
    - `YAMAGU32_YAMAGU34_...`: アンテナ1-2 相互相関
    - `YAMAGU66_YAMAGU66_...`: フェーズ合成後 自己相関
3.  **`..._geometricdelay.txt`**: 1秒ごとの幾何学的遅延の計算ログ。
4.  **`.png` プロット**: 相関強度、位相、スペクトル形状などを視覚化した画像。

## 遅延・レート補正の実装ノート（再現性重視）

この節は、実装時に実際に問題になった点を含めた「再現性チェックリスト」です。  
将来ソフトウェア FX 相関器を実装する際は、ここを仕様として固定してください。

### 1. まず固定すべき基準（ここが崩れると再現しない）
- **座標フレーム**: 入力 RA/Dec は **常に J2000** とみなす。
- **時刻基準**: `epoch` は UTC とし、内部は MJD/JD で扱う。
- **周波数基準**: 補正位相のキャリア周波数はアンテナごとのデータ帯域下端（`data_band_low_hz`）を使う。
- **符号規約**: 幾何遅延は `tau_g = -(r2-r1)·s / c`（`ant2-ant1`）。

過去に効いた修正:
- J2000 を of-date として扱っていたミスで rate 残差が増大。
- キャリア補正を `obsfreq` 固定でかけていたミスで残差が悪化。
- フレーム開始時刻評価により位相がずれるため、**フレーム中心時刻**評価へ統一。

### 2. 現行実装の補正モデル

幾何・clock・ユーザー補正を合成して、各フレームで `tau1/tau2` を計算します。

```text
tau_g(t)         : 幾何遅延（ant2-ant1）
tau_rel_no_clk(t)= tau_g(t) + gico3_correct + coarse + delay/fs + extra_rate*t
clock1(t)        = clock1_delay + clock1_rate*t
clock2(t)        = clock2_delay + clock2_rate*t

tau1(t) = (tau_rel_no_clk(t) - clock1(t)) - d_seek
tau2(t) = -clock2(t)
```

- `d_seek` は起動時の整数サンプル seek 分を秒換算したもの。
- 幾何 rate/accel は `tau_g(t)` の数値微分（中心差分）で得る。
- `--rate` は Hz で加算され、最終的に秒/秒へ換算して `tau_rel_no_clk(t)` に入る。

### 3. FX パイプラインで同一にすべき順序

`fringe` と `phased raw 出力` の両方で、以下を同じ順序で適用すること。

1. ビット展開（level/shuffle/sideband 正規化）
2. 整数遅延（time-domain サンプルシフト）
3. FFT
4. 小数遅延（周波数領域位相回転）
5. キャリア位相補正（`fr_lo1/fr_lo2`）
6. 共通帯域の bin 対応付け（band-align）
7. `X12 = X1 * conj(X2)` で相互相関、または重み付き加算して IFFT

重要:
- **fringe と synthesis で `tau1/tau2` 計算式は同一**であること。
- 差があると「相互相関は立つが phased raw が崩れる」状態になる。

### 4. なぜ相互相関スペクトルに 0 bin が出るか

この実装は「全帯域」ではなく **重複帯域のみ**を相関します。  
非重複 bin は意図的に未計算（0）です。これは FFT のゼロパディング由来ではなく、帯域マスク仕様です。

```text
band-align: shift N bins
band-overlap: M bins
```

この `M bins` 以外が 0 であれば正常です。

### 5. sideband/rotation の注意
- 入力 LSB は内部で USB に正規化してから処理。
- 出力 `.raw` は「両入力 LSB」のときのみ LSB に戻す。
- 大きな周波数オフセット（例: 343 MHz）は
  - 整数 bin シフト（band-align）と
  - 残差位相回転（rotation-fringestop residual）
  に分けて吸収する。

### 6. 再現性のための運用ルール

1. 必ず `cargo build --release` で実行する。  
2. 実行ログを保存し、次を毎回記録する。  
   `ra/dec`, `epoch`, `geom-delay/rate/accel`, `phase-freq`, `band-align`, `band-overlap`, `read-align delay`  
3. `frinZ` で `Res-Rate` を確認し、探索窓（通常 ±0.5 Hz）に収まることを確認する。  
4. 5秒積分など短時間では accel の寄与は小さい。先に delay/rate を詰める。  
5. 比較時は FFT 点数・sampling・帯域・sideband を固定する。

### 7. 将来のソフトウェア FX 相関器に向けた注意事項

- **座標・時刻・周波数の基準をコードで一元管理**する（分散実装しない）。
- 幾何モデル層と DSP 層を分離する。
  - 幾何: `tau(t), rate(t), accel(t)` を返すだけ
  - DSP: それを整数/小数遅延と位相回転へ写像
- 単体試験に「人工トーン」を必ず入れる。
  - 既知遅延/既知レートを注入し、回収できるか検証
- 出力検証を自動化する。
  - 自己相関ピーク位置
  - 相互相関の位相連続性
  - `Res-Rate` の上限閾値
- 高精度化するなら EOP（`dut1/xp/yp`）を導入する。
  - 現在は近似で十分なケースが多いが、長基線・長積分では効く。
