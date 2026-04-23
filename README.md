# Proton Transmitted Spectrum Calculation Through Cover Glass
## for Solar Cell Degradation Prediction Using NIEL

このREADMEは、`docs/report_proton_transmission.pdf` のページ構成（1/14〜14/14）に対応させて、
同じ流れで読めるように再構成したものです。  
図は `inputs/` を使ってコード実行した `outputs/` をそのまま使用しています。

## 実行条件と再現手順

### 入力ファイル（PDF対応）

- `inputs/proton_environment_LEO600km.csv`
  - Proton Energy (MeV)
  - Differential Fluence (cm-2 MeV-1)
- `inputs/range_borosilicate_glass.csv`
  - Energy (MeV)
  - Range_SRIM_NIST (um)
  - Range_Fitting (um)
  - Residual (%)

### 実行手順（Linux + uv）

```bash
cd /home/shirokawakita/Desktop/cursor/space_enviromental_analysis
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt
python src/proton_glass_transmission.py \
  inputs/proton_environment_LEO600km.csv \
  inputs/range_borosilicate_glass.csv
```

### 出力ファイル（outputs）

- `outputs/fig3_range_fit.csv`
- `outputs/fig3_range_fit.png`
- `outputs/transmitted_fluence.csv`
- `outputs/transmitted_fluence.png`

---

## PDF 1/14 対応

**タイトルページ**

- Proton Transmitted Spectrum Calculation Through Cover Glass
- for Solar Cell Degradation Prediction Using NIEL

本READMEは上記タイトルの技術資料内容を、実装付きで再現しています。

## PDF 2/14 対応

**DDD評価式（導入）**

- `Dp = Φp(E) · Sp(E) [MeV/g]`
- `Dp = ∫ (dΦp(E)/dE) · Sp(E) dE`
- `Dd_total = Dp + De,eff(1.0) / Rep`

本コードの役割は、この評価で使う陽子項の入力として、透過後スペクトル `f(ε)` を生成することです。

## PDF 3/14 対応

**飛程近似式（式1）**

- `R(E) = A · E^a + B · E^b [um]`
- `A = 202.0955`
- `a = 1.5279`
- `B = 14.7058`
- `b = 1.0166`

この飛程式を `src/proton_glass_transmission.py` 内で使用し、逆関数計算で必要入射エネルギーを求めます。

## PDF 4/14 対応

**図面・説明ページ（資料構成上の中間ページ）**

PDFではこのページは式のつなぎとして使われています。  
実装側では、次ページ以降の透過条件式と角度積分準備に対応します。

## PDF 5/14 対応

**透過後エネルギー関係（式2）**

- `R(ε) = R(E) - t / cosθ`

厚み条件（PDF記載）:

- 0.05 mm -> 50 um -> 最小透過エネルギー 約0.372 MeV
- 0.10 mm -> 100 um -> 最小透過エネルギー 約0.595 MeV
- 0.15 mm -> 150 um -> 最小透過エネルギー 約0.781 MeV
- 0.20 mm -> 200 um -> 最小透過エネルギー 約0.947 MeV

## PDF 6/14 対応

**図面・説明ページ（資料構成上の中間ページ）**

実装では、透過後エネルギーグリッドと角度離散化の設定に対応します。

## PDF 7/14 対応

**スペクトル変換の導出（式3, 式4）**

- `f(ε) dε = g(E) dE`
- `f(ε) = g(E) · |dE/dε|`
- `R'(ε) = R'(E) · dE/dε`
- `dE/dε = R'(ε)/R'(E)`

ここで `g(E)` は入力（4π）スペクトルです。

## PDF 8/14 対応

**最終積分式（式5）**

- `f(ε) = (1/2) ∫ g[E(θ)] · R'(ε)/R'(E(θ)) · sinθ dθ`

補足:

- `1/2` は 4π -> 2π 変換係数
- 角度重みは `sinθ`
- 実装では `θ = 0` から `π/2` 直前までを離散化積分

## PDF 9/14 対応

**図面・説明ページ（資料構成上の中間ページ）**

実装では、厚みごとの透過後スペクトル計算ループに対応します。

## PDF 10/14 対応

**スペクトル結果ページ（図）**

コード出力における対応図:

![Transmitted Fluence](outputs/transmitted_fluence.png)

読み方:

- 黒丸: Uncovered (4π)
- 実線/破線: 50/100/150/200 um の透過後スペクトル (2π)
- 厚み増加で低エネルギー側が抑制される

## PDF 11/14 対応

**代表値（積分フルエンス）**

- 50 um (0.05 mm): `4.41 x 10^10 cm^-2` (1.00, reference)
- 100 um (0.10 mm): `3.07 x 10^10 cm^-2` (0.70)
- 150 um (0.15 mm): `2.42 x 10^10 cm^-2` (0.55)
- 200 um (0.20 mm): `2.02 x 10^10 cm^-2` (0.46)

`outputs/transmitted_fluence.csv` から同傾向を確認できます。

## PDF 12/14 対応

**計算条件一覧**

- `GLASS_THICKNESSES_MM = [0.05, 0.10, 0.15, 0.20]`
- `N_ENERGY_POINTS = 200`
- `N_THETA = 200`
- Range fitting params: `A, a, B, b`
- 実行コマンド: `python src/proton_glass_transmission.py <environment_csv> <range_csv>`

## PDF 13/14 対応

**参考文献**

PDFの参考文献ページに準拠し、NIEL/DDD関連文献を参照しています。

## PDF 14/14 対応

**終端ページ**

資料末尾に対応します。

---

## Fig.3相当（Range fit）

PDF内の飛程式検証に対応する出力図:

![Fig3 Range Fit](outputs/fig3_range_fit.png)

## `src/proton_glass_transmission.py` の詳細解説

### 処理の全体フロー

1. 2つのCSV（環境スペクトル、飛程参照）を読み込み
2. 環境スペクトル `g(E)` を対数補間関数化
3. 飛程フィット妥当性（Fig.3）を `csv/png` 出力
4. 透過後エネルギー軸 `ε` を作成
5. 各ガラス厚に対して式5を離散角度積分で評価
6. 透過後スペクトルを `csv/png` 出力

### 主要定数・設定値

- `GLASS_THICKNESSES_MM`: 解析対象厚み（mm）
- `N_ENERGY_POINTS`: 透過後エネルギー点数（精度と速度のトレードオフ）
- `N_THETA`: 角度積分分割数（精度と速度のトレードオフ）
- `EPS_MIN`, `EPS_MAX`: 透過後エネルギー範囲
- `RANGE_A`, `RANGE_a`, `RANGE_B`, `RANGE_b`: 飛程近似パラメータ

### 関数ごとの役割

- `range_func(E)`
  - 式1 `R(E)` を返すコア関数
  - スカラー/配列のどちらにも対応
- `dR_dE(E)`
  - `R(E)` の導関数 `R'(E)` を返す
  - ヤコビアン `dE/dε` の算出に使用
- `find_energy_from_range(R_target, E_max=2e4)`
  - 逆関数 `R(E)=R_target` を `brentq` で解く
  - 式2から入射エネルギーを復元する要所

- `read_environment_csv(filepath)`
  - `inputs/proton_environment_LEO600km.csv` 用の読込
  - 正値のみ採用し、エネルギー順にソート
- `read_range_reference_csv(filepath)`
  - `inputs/range_borosilicate_glass.csv` 用の読込
  - 1列目Energy・2列目Rangeを主に利用

- `create_fluence_interpolator(E_env, F_env)`
  - `log10(E)-log10(F)` 線形補間で `g(E)` を生成
  - 範囲外は0扱い（外挿を避ける安全設計）

- `compute_transmitted_spectrum(t_glass_um, eps_array, diff_fluence_func, N_theta, verbose)`
  - 本コードの中心処理
  - 各 `ε` について各 `θ` を走査し、次を積算:
    - `path_length = t/cosθ`
    - `R_required = R(ε) + path_length`
    - `E_inc = R^-1(R_required)`
    - `jacobian = R'(ε)/R'(E_inc)`
    - `weight = sinθ * dθ`
    - 寄与 `g(E_inc) * jacobian * weight`
  - 最後に `0.5` を掛けて 4π -> 2π 変換

- `write_fig3_csv(...)`, `create_fig3_plot(...)`
  - 飛程参照値とフィット値を比較し、残差を含むFig.3相当を出力
- `write_output_csv(...)`, `create_plot(...)`
  - 透過後スペクトル結果を表と図で出力

- `main()`
  - 引数検証、`outputs/` 作成、読込、計算、出力までのオーケストレーション

### 数値計算上の注意点

- `θ = π/2` では `cosθ -> 0` で発散するため、実装では `π/2 - 0.005` まで積分
- `E_inc <= ε` や `g(E)<=0` は物理的寄与なしとしてスキップ
- 積分フルエンスは `np.trapezoid` で計算し、厚み依存の減衰を定量化

### 計算量と性能感

- 計算量は概ね `O(厚み数 * N_ENERGY_POINTS * N_THETA * root_solve_cost)`
- デフォルト（厚み4, エネルギー200, 角度200）で数秒程度
- 精度を上げるときは `N_THETA` と `N_ENERGY_POINTS` を段階的に増やすのが安全

### 拡張ポイント

- ガラス厚追加: `GLASS_THICKNESSES_MM` を変更
- 分解能調整: `N_ENERGY_POINTS`, `N_THETA` を変更
- 材料変更: `RANGE_A/a/B/b` を対象材料のフィット値に差し替え
- 入力環境変更: `inputs/proton_environment_LEO600km.csv` を別軌道データへ差し替え

### 関数名・PDFページ・式番号の対応

- `range_func(E)`
  - PDF対応: 3/14
  - 式対応: 式1（飛程近似）
- `dR_dE(E)`
  - PDF対応: 7/14, 8/14
  - 式対応: 式4, 式5 内の `R'(E)` / `R'(ε)`
- `find_energy_from_range(R_target, E_max)`
  - PDF対応: 5/14
  - 式対応: 式2（`R(ε) = R(E) - t/cosθ` の逆算）
- `read_environment_csv(filepath)`
  - PDF対応: 10/14（Uncoveredスペクトル入力）
  - 式対応: `g(E)` 入力データ
- `read_range_reference_csv(filepath)`
  - PDF対応: 3/14, Fig.3検証
  - 式対応: 式1の妥当性確認用参照データ
- `create_fluence_interpolator(E_env, F_env)`
  - PDF対応: 7/14, 8/14
  - 式対応: `g(E)` の連続評価
- `compute_transmitted_spectrum(...)`
  - PDF対応: 5/14, 7/14, 8/14, 9/14
  - 式対応: 式2 -> 式4 -> 式5 の主計算
- `write_fig3_csv(...)`, `create_fig3_plot(...)`
  - PDF対応: 3/14（飛程式）とFig.3
  - 式対応: 式1の参照値比較（残差）
- `write_output_csv(...)`, `create_plot(...)`
  - PDF対応: 10/14, 11/14
  - 式対応: 式5の計算結果出力
- `main()`
  - PDF対応: 全ページ統合
  - 式対応: 全式の実行順制御

## 依存モジュール

- `numpy`
- `scipy`
- `matplotlib`

## 関連資料

- 元資料: `docs/report_proton_transmission.pdf`
- 実装詳細: `docs/methodology_proton_transmission.md`
