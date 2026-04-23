#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==========================================================================
 ホウケイ酸ガラス（カバーガラス）透過後の陽子線微分フルエンス計算
 Proton Differential Fluence After Borosilicate Glass Transmission
==========================================================================

 概要:
   軌道上の捕捉陽子線環境（4π空間微分フルエンス）を入力として、
   ホウケイ酸ガラス透過後の微分フルエンス（2π空間）を計算する。

 理論:
   - 飛程: R(E) = A*E^a + B*E^b  (式46)
   - 透過後エネルギー: R(ε) = R(E) - t/cosθ  (式50)
   - 透過後スペクトル: f(ε) = Σ_θ g(E) * (dE/dε) * sinθ * Δθ  (式52)

入力:
  1) 環境CSV
     A列: Proton Energy [MeV]
     B列: Differential Fluence [cm^-2 MeV^-1]（4π空間）
  2) 飛程検証CSV（Fig.3 用）
     A列: Energy [MeV]
     B列: Range_SRIM_NIST [μm]

出力:
  - CSV: 透過後微分フルエンス（各ガラス厚, 図10相当）
  - PNG: スペクトル比較プロット（図10相当）
  - CSV: 飛程フィット比較（図3相当）
  - PNG: 飛程フィット比較（図3相当）

 参考文献:
   [1] S.R. Messenger et al., Prog. Photovoltaics 5 (1997) 407-413.
   [2] S.R. Messenger et al., IEEE TNS 44 (1997) 2169-2173.

使用方法:
  python proton_glass_transmission.py <environment_csv> <range_csv>

 必要なライブラリ:
   numpy, scipy, matplotlib
==========================================================================
"""

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')  # GUIバックエンド不要
import matplotlib.pyplot as plt
import csv
import sys
import os
import time

# =====================================================================
# 設定パラメータ（必要に応じて変更）
# =====================================================================

# ガラス厚 [mm] のリスト
GLASS_THICKNESSES_MM = [0.05, 0.10, 0.15, 0.20]

# 透過後エネルギーグリッド数（多いほど精度向上、計算時間増加）
N_ENERGY_POINTS = 200

# 角度分割数（多いほど精度向上、計算時間増加）
N_THETA = 200

# 透過後エネルギー範囲 [MeV]
EPS_MIN = 1e-4
EPS_MAX = 1e3

# 出力ディレクトリとファイル名
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs")
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "transmitted_fluence.csv")
OUTPUT_PNG = os.path.join(OUTPUT_DIR, "transmitted_fluence.png")
OUTPUT_FIG3_CSV = os.path.join(OUTPUT_DIR, "fig3_range_fit.csv")
OUTPUT_FIG3_PNG = os.path.join(OUTPUT_DIR, "fig3_range_fit.png")

# =====================================================================
# ホウケイ酸ガラス中の陽子線飛程パラメータ
# R(E) = A * E^a + B * E^b   [μm]  (E in MeV)
# SRIM/NIST データからのフィッティング
# =====================================================================
RANGE_A = 202.0955
RANGE_a = 1.5279
RANGE_B = 14.7058
RANGE_b = 1.0166


def range_func(E):
    """
    ホウケイ酸ガラス中の陽子線飛程 R(E) [μm]

    Parameters:
        E: エネルギー [MeV] (scalar or array)

    Returns:
        R: 飛程 [μm]
    """
    return RANGE_A * np.power(E, RANGE_a) + RANGE_B * np.power(E, RANGE_b)


def dR_dE(E):
    """
    飛程の微分 dR/dE [μm/MeV]

    Parameters:
        E: エネルギー [MeV]

    Returns:
        dR/dE
    """
    return (RANGE_a * RANGE_A * np.power(E, RANGE_a - 1)
            + RANGE_b * RANGE_B * np.power(E, RANGE_b - 1))


def find_energy_from_range(R_target, E_max=2e4):
    """
    飛程 R_target [μm] に対応するエネルギー E [MeV] を数値的に求める（逆関数）

    Parameters:
        R_target: 目標飛程 [μm]
        E_max: 探索上限エネルギー [MeV]

    Returns:
        E: エネルギー [MeV]、見つからない場合は None
    """
    if R_target <= 0:
        return None
    try:
        return brentq(lambda E: range_func(E) - R_target, 1e-8, E_max, xtol=1e-12)
    except (ValueError, RuntimeError):
        return None


def read_environment_csv(filepath):
    """
    環境データCSVを読み込む

    Parameters:
        filepath: CSVファイルパス
            A列: Proton Energy [MeV]
            B列: Differential Fluence [cm^-2 MeV^-1]

    Returns:
        E_env: エネルギー配列 [MeV]
        F_env: 微分フルエンス配列 [cm^-2 MeV^-1]
    """
    E_list = []
    F_list = []

    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)

        for row in reader:
            if len(row) < 2:
                continue
            try:
                e_val = float(row[0])
                f_val = float(row[1])
                if e_val > 0 and f_val > 0:
                    E_list.append(e_val)
                    F_list.append(f_val)
            except (ValueError, IndexError):
                continue  # ヘッダ行や不正データはスキップ

    if len(E_list) < 2:
        raise ValueError("CSVファイルに有効なデータが2点以上必要です。")

    E_env = np.array(E_list)
    F_env = np.array(F_list)

    # エネルギー順にソート
    sort_idx = np.argsort(E_env)
    E_env = E_env[sort_idx]
    F_env = F_env[sort_idx]

    return E_env, F_env


def read_range_reference_csv(filepath):
    """
    飛程検証用CSVを読み込む（Fig.3 用）

    Parameters:
        filepath: CSVファイルパス
            A列: Energy [MeV]
            B列: Range_SRIM_NIST [μm]

    Returns:
        E_ref: エネルギー配列 [MeV]
        R_ref: 参照飛程配列 [μm]
    """
    E_list = []
    R_list = []

    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)

        for row in reader:
            if len(row) < 2:
                continue
            try:
                e_val = float(row[0])
                r_val = float(row[1])
                if e_val > 0 and r_val > 0:
                    E_list.append(e_val)
                    R_list.append(r_val)
            except (ValueError, IndexError):
                continue

    if len(E_list) < 2:
        raise ValueError("飛程CSVファイルに有効なデータが2点以上必要です。")

    E_ref = np.array(E_list)
    R_ref = np.array(R_list)

    sort_idx = np.argsort(E_ref)
    E_ref = E_ref[sort_idx]
    R_ref = R_ref[sort_idx]

    return E_ref, R_ref


def create_fluence_interpolator(E_env, F_env):
    """
    微分フルエンスの対数補間関数を作成

    Parameters:
        E_env: エネルギー配列 [MeV]
        F_env: 微分フルエンス配列 [cm^-2 MeV^-1]

    Returns:
        interp_func: 補間関数 E -> F（範囲外は0を返す）
    """
    log_interp = interp1d(
        np.log10(E_env), np.log10(F_env),
        kind='linear', fill_value=-np.inf, bounds_error=False
    )

    def interp_func(E):
        E = np.atleast_1d(np.float64(E))
        result = np.zeros_like(E)
        mask = (E >= E_env.min()) & (E <= E_env.max())
        if mask.any():
            log_f = log_interp(np.log10(E[mask]))
            result[mask] = np.power(10.0, log_f)
        return result

    return interp_func


def compute_transmitted_spectrum(t_glass_um, eps_array, diff_fluence_func,
                                  N_theta=200, verbose=True):
    """
    ガラス透過後の微分フルエンスを計算（式52-53に基づく）

    Parameters:
        t_glass_um: ガラス厚 [μm]
        eps_array: 透過後エネルギー配列 [MeV]
        diff_fluence_func: 入射微分フルエンス関数 g(E) [cm^-2 MeV^-1]（4π空間）
        N_theta: 角度分割数
        verbose: 進捗表示

    Returns:
        f_eps: 透過後微分フルエンス [cm^-2 MeV^-1]（2π空間）

    Note:
        - 入力は4π空間のフルエンス、出力は2π空間（裏面半無限遮蔽）
        - 等方入射を仮定
        - 式(52): f(ε) = (1/2) × Σ_k g[E(θ_k)] × R'(ε)/R'(E) × sinθ_k × Δθ
        - 1/2 は 4π→2π（半球）変換係数
    """
    theta_max = np.pi / 2 - 0.005  # π/2 直前まで
    theta_arr = np.linspace(0, theta_max, N_theta)
    d_theta = theta_arr[1] - theta_arr[0]

    f_eps = np.zeros(len(eps_array))
    n_total = len(eps_array)
    t_start = time.time()

    for i, eps in enumerate(eps_array):
        if eps < 1e-8:
            continue

        R_eps = range_func(eps)
        dRdE_at_eps = dR_dE(eps)

        total = 0.0

        for theta in theta_arr:
            # 遮蔽材中のパス長
            path_length = t_glass_um / np.cos(theta)

            # 必要な飛程
            R_required = R_eps + path_length

            # 入射エネルギーを逆算
            E_inc = find_energy_from_range(R_required)
            if E_inc is None or E_inc <= eps:
                continue

            # 入射スペクトル g(E)（4π空間の値をそのまま使用）
            g_E = diff_fluence_func(np.array([E_inc]))[0]
            if g_E <= 0:
                continue

            # ヤコビアン dE/dε = R'(ε) / R'(E)
            # 導出: R(ε) = R(E) - t/cosθ → R'(ε) = R'(E)·dE/dε
            #       ∴ dE/dε = R'(ε) / R'(E)  (<1, スペクトル圧縮)
            dRdE_at_E = dR_dE(E_inc)
            jacobian = dRdE_at_eps / dRdE_at_E

            # 立体角の重み
            weight = np.sin(theta) * d_theta

            total += g_E * jacobian * weight

        # 4π→2π 変換: 等方入射の半球積分は全球の1/2
        f_eps[i] = total * 0.5

        # 進捗表示
        if verbose and ((i + 1) % 50 == 0 or i == n_total - 1):
            elapsed = time.time() - t_start
            print(f"    t={t_glass_um:6.0f} μm: {i+1:4d}/{n_total} "
                  f"({(i+1)/n_total*100:5.1f}%) "
                  f"[{elapsed:.1f}s]", flush=True)

    return f_eps


def write_fig3_csv(filepath, E_ref, R_ref, R_fit, residual_pct):
    """
    Fig.3 相当の飛程フィット比較をCSV出力
    """
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            "Energy (MeV)",
            "Range_SRIM_NIST (um)",
            "Range_Fitting (um)",
            "Residual (%)"
        ])
        for i in range(len(E_ref)):
            writer.writerow([
                f"{E_ref[i]:.6E}",
                f"{R_ref[i]:.6E}",
                f"{R_fit[i]:.6E}",
                f"{residual_pct[i]:.4f}"
            ])


def write_output_csv(filepath, eps_array, uncov_2pi, results_dict, thicknesses_um):
    """
    計算結果をCSVファイルに出力

    Parameters:
        filepath: 出力ファイルパス
        eps_array: エネルギー配列 [MeV]
        uncov_2pi: 遮蔽なし2π空間微分フルエンス
        results_dict: {厚さ[μm]: 透過後微分フルエンス配列}
        thicknesses_um: ガラス厚リスト [μm]
    """
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        # ヘッダ
        header = ["Energy (MeV)", "Uncovered 2pi (cm-2 MeV-1)"]
        for t_um in thicknesses_um:
            header.append(f"{t_um:.0f}um ({t_um/1000:.2f}mm) 2pi (cm-2 MeV-1)")
        writer.writerow(header)

        # データ
        for i in range(len(eps_array)):
            row = [f"{eps_array[i]:.6E}", f"{uncov_2pi[i]:.6E}"]
            for t_um in thicknesses_um:
                row.append(f"{results_dict[t_um][i]:.6E}")
            writer.writerow(row)


def create_fig3_plot(filepath, E_ref, R_ref, R_fit, residual_pct):
    """
    Fig.3 相当: 飛程データとフィット、および残差を可視化
    """
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 8), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]}
    )

    ax1.plot(E_ref, R_ref, 'ko', markersize=4, label='SRIM/NIST data')
    ax1.plot(E_ref, R_fit, '-', color='black', linewidth=1.4, label='Fitting')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel('Range, μm', fontsize=12)
    ax1.grid(True, which='major', alpha=0.3)
    ax1.grid(True, which='minor', alpha=0.1)
    ax1.legend(loc='upper left', fontsize=10, frameon=True, fancybox=False)

    ax2.plot(E_ref, residual_pct, 's-', color='black', markersize=3, linewidth=1.0)
    ax2.axhline(0.0, color='gray', linestyle='--', linewidth=1.0)
    ax2.set_xscale('log')
    ax2.set_xlabel('Proton Energy, MeV', fontsize=12)
    ax2.set_ylabel('Residual, %', fontsize=12)
    ax2.grid(True, which='major', alpha=0.3)
    ax2.grid(True, which='minor', alpha=0.1)

    plt.tight_layout()
    plt.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Fig.3プロット保存: {filepath}")


def create_plot(filepath, eps_array, E_env, F_env, uncov_2pi,
                results_dict, thicknesses_um):
    """
    スペクトル比較プロット（図10相当）を作成

    Parameters:
        filepath: 出力PNGファイルパス
        eps_array: エネルギー配列 [MeV]
        E_env: 入力環境データのエネルギー [MeV]
        F_env: 入力環境データの微分フルエンス [cm^-2 MeV^-1]
        uncov_2pi: 遮蔽なし2πフルエンス
        results_dict: 計算結果
        thicknesses_um: ガラス厚リスト [μm]
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    # 入力環境データ（4π、プロット点）
    ax.plot(E_env, F_env, 'ko', markersize=4, label='Uncovered (4π)', zorder=10)

    # 線スタイルの定義
    line_styles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1, 1, 1))]
    colors = ['black', '#333333', '#555555', '#777777', '#999999']

    # 各ガラス厚の透過スペクトル
    for idx, t_um in enumerate(thicknesses_um):
        f_trans = results_dict[t_um]
        mask = f_trans > 0
        if mask.any():
            style = line_styles[idx % len(line_styles)]
            color = colors[idx % len(colors)]
            ax.plot(
                eps_array[mask], f_trans[mask],
                linestyle=style, color=color, linewidth=1.2,
                label=f'{t_um:.0f}μm ({t_um/1000:.2f}mm) (2π)'
            )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Proton Energy, MeV', fontsize=13)
    ax.set_ylabel('Differential Fluence, cm$^{-2}$MeV$^{-1}$', fontsize=13)
    ax.set_title('Proton Transmitted Spectrum Through Cover Glass', fontsize=13, pad=10)
    ax.set_xlim(1e-4, 1e3)

    # Y軸範囲を自動調整
    all_max = max(F_env.max(), max(results_dict[t].max() for t in thicknesses_um))
    all_min_list = [f_trans[f_trans > 0].min() for t, f_trans in results_dict.items()
                    if (f_trans > 0).any()]
    all_min = min(all_min_list) if all_min_list else 1e0
    ax.set_ylim(all_min * 0.1, all_max * 5)

    ax.legend(loc='upper right', fontsize=9, frameon=True, fancybox=False,
              edgecolor='black', framealpha=0.9)
    ax.grid(True, which='major', alpha=0.3)
    ax.grid(True, which='minor', alpha=0.1)

    # 軌道情報テキスト（カスタマイズ可能）
    ax.text(0.72, 0.55, 'Coverglass thickness\n(Borosilicate glass)',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

    ax.tick_params(which='both', direction='in', top=True, right=True)

    plt.tight_layout()
    plt.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"\n  プロット保存: {filepath}")


# =====================================================================
# メイン処理
# =====================================================================
def main():
    print("=" * 70)
    print("  ホウケイ酸ガラス透過後の陽子線微分フルエンス計算")
    print("  Proton Transmitted Spectrum Through Borosilicate Glass")
    print("=" * 70)

    # -----------------------------------------------------------------
    # 引数チェック
    # -----------------------------------------------------------------
    if len(sys.argv) < 3:
        print("\n使用方法:")
        print(f"  python {sys.argv[0]} <environment_csv> <range_csv>")
        print("\n環境CSV形式:")
        print("  A列: Proton Energy [MeV]")
        print("  B列: Differential Fluence [cm^-2 MeV^-1] (4π空間)")
        print("\n飛程CSV形式:")
        print("  A列: Energy [MeV]")
        print("  B列: Range_SRIM_NIST [μm]")
        print("\nヘッダ行は自動的にスキップされます。")
        sys.exit(1)

    input_csv = sys.argv[1]
    range_csv = sys.argv[2]

    if not os.path.exists(input_csv):
        print(f"\nエラー: ファイル '{input_csv}' が見つかりません。")
        sys.exit(1)
    if not os.path.exists(range_csv):
        print(f"\nエラー: ファイル '{range_csv}' が見つかりません。")
        sys.exit(1)

    # 出力先ディレクトリ作成
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\n出力ディレクトリ: {OUTPUT_DIR}")

    # -----------------------------------------------------------------
    # 環境データ読み込み
    # -----------------------------------------------------------------
    print(f"\n入力ファイル: {input_csv}")
    E_env, F_env = read_environment_csv(input_csv)
    print(f"  データ点数: {len(E_env)}")
    print(f"  エネルギー範囲: {E_env.min():.4E} ~ {E_env.max():.4E} MeV")
    print(f"  フルエンス範囲: {F_env.min():.4E} ~ {F_env.max():.4E} cm^-2 MeV^-1")

    print(f"\n飛程ファイル: {range_csv}")
    E_ref, R_ref = read_range_reference_csv(range_csv)
    print(f"  データ点数: {len(E_ref)}")
    print(f"  エネルギー範囲: {E_ref.min():.4E} ~ {E_ref.max():.4E} MeV")
    print(f"  飛程範囲: {R_ref.min():.4E} ~ {R_ref.max():.4E} μm")

    # 補間関数の作成
    diff_fluence_func = create_fluence_interpolator(E_env, F_env)

    # -----------------------------------------------------------------
    # ガラス厚の設定
    # -----------------------------------------------------------------
    thicknesses_um = [t * 1000.0 for t in GLASS_THICKNESSES_MM]  # mm → μm
    print(f"\nガラス厚: {GLASS_THICKNESSES_MM} mm")
    print(f"         = {thicknesses_um} μm")

    # -----------------------------------------------------------------
    # 飛程パラメータの表示
    # -----------------------------------------------------------------
    print(f"\n飛程フィッティング: R(E) = {RANGE_A:.4f}*E^{RANGE_a:.4f} "
          f"+ {RANGE_B:.4f}*E^{RANGE_b:.4f} [μm]")

    for t_um in thicknesses_um:
        E_min = find_energy_from_range(t_um)
        if E_min:
            print(f"  t={t_um:.0f}μm: 垂直入射最小透過エネルギー = {E_min:.4f} MeV")

    # -----------------------------------------------------------------
    # Fig.3相当: 飛程フィット比較の生成
    # -----------------------------------------------------------------
    R_fit = range_func(E_ref)
    residual_pct = (R_fit - R_ref) / R_ref * 100.0
    write_fig3_csv(OUTPUT_FIG3_CSV, E_ref, R_ref, R_fit, residual_pct)
    print(f"  Fig.3 CSV保存: {OUTPUT_FIG3_CSV}")
    create_fig3_plot(OUTPUT_FIG3_PNG, E_ref, R_ref, R_fit, residual_pct)

    # -----------------------------------------------------------------
    # 透過後エネルギーグリッドの生成
    # -----------------------------------------------------------------
    eps_array = np.logspace(np.log10(EPS_MIN), np.log10(EPS_MAX), N_ENERGY_POINTS)
    print(f"\n透過後エネルギーグリッド: {N_ENERGY_POINTS} 点")
    print(f"  範囲: {EPS_MIN:.4E} ~ {EPS_MAX:.4E} MeV")
    print(f"角度分割数: {N_theta}")

    # -----------------------------------------------------------------
    # 遮蔽なし（2π空間 = 4πの半分）
    # -----------------------------------------------------------------
    uncov_2pi = diff_fluence_func(eps_array) / 2.0

    # -----------------------------------------------------------------
    # 各ガラス厚での計算
    # -----------------------------------------------------------------
    results_dict = {}
    t_total_start = time.time()

    for t_um in thicknesses_um:
        print(f"\n--- ガラス厚 {t_um:.0f} μm ({t_um/1000:.2f} mm) ---")
        f_trans = compute_transmitted_spectrum(
            t_um, eps_array, diff_fluence_func,
            N_theta=N_THETA, verbose=True
        )
        results_dict[t_um] = f_trans

        # 積分フルエンスの表示
        mask = f_trans > 0
        if mask.any():
            total = np.trapezoid(f_trans, eps_array)
            print(f"  積分フルエンス: {total:.3E} cm^-2")

    elapsed_total = time.time() - t_total_start
    print(f"\n全計算時間: {elapsed_total:.1f} 秒")

    # -----------------------------------------------------------------
    # CSV出力
    # -----------------------------------------------------------------
    print(f"\n--- 出力 ---")
    write_output_csv(OUTPUT_CSV, eps_array, uncov_2pi, results_dict, thicknesses_um)
    print(f"  CSV保存: {OUTPUT_CSV}")

    # -----------------------------------------------------------------
    # PNG出力
    # -----------------------------------------------------------------
    create_plot(OUTPUT_PNG, eps_array, E_env, F_env, uncov_2pi,
                results_dict, thicknesses_um)

    print(f"\n{'=' * 70}")
    print("  完了")
    print(f"{'=' * 70}")


# グローバル変数の参照（compute_transmitted_spectrum 内で使用）
N_theta = N_THETA

if __name__ == "__main__":
    main()
