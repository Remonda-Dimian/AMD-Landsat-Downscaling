#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AMD Validation — RGB only.
Implements Grad.RI (7), Freq.RI (8), SSIM.RI (9), Weighted.RI (10).
"""

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from scipy import fft
import cv2
from skimage.metrics import structural_similarity as ssim_sk
import warnings
warnings.filterwarnings('ignore')

BAND_MAP = {
    'B2': {'l_band': 1, 's2_band': 1, 'name': 'Blue'},
    'B3': {'l_band': 2, 's2_band': 2, 'name': 'Green'},
    'B4': {'l_band': 3, 's2_band': 3, 'name': 'Red'},
}

def read_and_align(filepath, band_idx, target_shape, target_transform, target_crs):
    with rasterio.open(filepath) as src:
        dst = np.empty(target_shape, dtype=np.float32)
        reproject(
            source=rasterio.band(src, band_idx), destination=dst,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=target_transform, dst_crs=target_crs,
            resampling=Resampling.bilinear
        )
    return dst

def compute_grad_ri(d, r):
    gx_d = cv2.Sobel(d, cv2.CV_32F, 1, 0, ksize=3)
    gy_d = cv2.Sobel(d, cv2.CV_32F, 0, 1, ksize=3)
    gx_r = cv2.Sobel(r, cv2.CV_32F, 1, 0, ksize=3)
    gy_r = cv2.Sobel(r, cv2.CV_32F, 0, 1, ksize=3)
    
    sum_d = np.nansum(np.sqrt(gx_d**2 + gy_d**2))
    sum_r = np.nansum(np.sqrt(gx_r**2 + gy_r**2))
    
    return sum_d / sum_r if sum_r != 0 else np.nan

def compute_freq_ri(d, r):
    F_d = fft.fftshift(fft.fft2(d))
    F_r = fft.fftshift(fft.fft2(r))
    
    rows, cols = d.shape
    crow, ccol = rows // 2, cols // 2
    u = (np.arange(cols) - ccol) / (cols / 2)
    v = (np.arange(rows) - crow) / (rows / 2)
    U, V = np.meshgrid(u, v)
    hf = np.sqrt(U**2 + V**2) > 0.15
    
    sum_d = np.nansum(np.abs(F_d)[hf])
    sum_r = np.nansum(np.abs(F_r)[hf])
    
    return sum_d / sum_r if sum_r != 0 else np.nan

def compute_ssim_ri(d, r):
    d = np.clip(d, 0, 1)
    r = np.clip(r, 0, 1)
    return ssim_sk(d, r, data_range=1.0, win_size=7, channel_axis=None)

def compute_weighted_ri(g, f, s):
    return 0.4 * g + 0.3 * f + 0.3 * s

def validate_band(landsat_path, s2_path, band_key, cloud_mask_path=None):
    cfg = BAND_MAP[band_key]
    
    with rasterio.open(s2_path) as src:
        s2_t = src.transform
        s2_shape = (src.height, src.width)
        s2_crs = src.crs
        s2_data = src.read(cfg['s2_band'])
    
    l_data = read_and_align(landsat_path, cfg['l_band'], s2_shape, s2_t, s2_crs)
    
    if cloud_mask_path:
        with rasterio.open(cloud_mask_path) as src:
            cmask = src.read(1)
            if cmask.shape != s2_shape:
                cmask = cv2.resize(cmask, (s2_shape[1], s2_shape[0]), interpolation=cv2.INTER_NEAREST)
        valid = (cmask == 1) & np.isfinite(s2_data) & np.isfinite(l_data) & (s2_data > 0) & (l_data > 0)
    else:
        valid = np.isfinite(s2_data) & np.isfinite(l_data) & (s2_data > 0) & (l_data > 0)
    
    d = np.where(valid, l_data, np.nan)
    r = np.where(valid, s2_data, np.nan)
    
    return {
        'band': band_key,
        'name': cfg['name'],
        'Grad_RI': float(compute_grad_ri(d, r)),
        'Freq_RI': float(compute_freq_ri(d, r)),
        'SSIM_RI': float(compute_ssim_ri(d, r)),
        'Weighted_RI': float(compute_weighted_ri(
            compute_grad_ri(d, r), compute_freq_ri(d, r), compute_ssim_ri(d, r)
        ))
    }

def main():
    landsat = '/path/to/amd_landsat_10m_rgb.tif'
    s2 = '/path/to/sentinel2_10m_rgb.tif'
    cmask = '/path/to/cloud_mask.tif'
    
    print("AMD Validation — RGB")
    print("=" * 60)
    
    results = []
    for b in BAND_MAP.keys():
        print(f"\n>>> {BAND_MAP[b]['name']} ({b})...")
        try:
            res = validate_band(landsat, s2, b, cmask)
            results.append(res)
            for k in ['Grad_RI', 'Freq_RI', 'SSIM_RI', 'Weighted_RI']:
                print(f"    {k:<12}: {res[k]:.4f}")
        except Exception as e:
            print(f"    ERROR: {e}")
    
    print("\n" + "=" * 60)
    print(f"{'Band':<8} {'Name':<10} {'Grad.RI':<10} {'Freq.RI':<10} {'SSIM.RI':<10} {'W.RI':<10}")
    print("-" * 60)
    for r in results:
        print(f"{r['band']:<8} {r['name']:<10} {r['Grad_RI']:<10.4f} "
              f"{r['Freq_RI']:<10.4f} {r['SSIM_RI']:<10.4f} {r['Weighted_RI']:<10.4f}")
    print("=" * 60)

if __name__ == "__main__":
    main()
